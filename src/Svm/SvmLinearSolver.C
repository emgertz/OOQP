#include "sparseutils.h"
#include "SvmLinearSolver.h"
#include "OoqpBlas.h"
#include "SimpleVector.h"
#include <math.h>
#include <string.h>

#include<vector>
#include<algorithm>

extern int gOoqpPrintLevel;
int precon_count = 0;

// A factor by which to scale the svmtol.
double SvmPrecond::  mScalesvmtol      =  100;

extern "C" {
  void dspr_( char *, int *, double *, double *, int *, double * );

  // void GetTime( double  * utime, double * stime );

  void dpptrf_( char * uplo, int * n,
                double * A, int * info );

  void dpptrs_( char * uplo, int * n, int * nrhs,
                double * A,
                double * b, int * ldb, int * info ) ;

  void daxpy_ ( int * n, double * alpha, double * x,
                int * incx, double * y, int * incy );
  double ddot_( int * n, double * dx, int * incx, double * dy, int * incy);
  double dnrm2_ (int * n, double * x, int * incx );
}


//////////////////////////////////////////////////////////////////////
SvmPrecond::SvmPrecond( int m, int n, int usingDirectSolve ) :
  mP(0), mPfactored(0), nobs(m), hdim(n), mNsv(0), mGammaReduced(0),
  mUsingDirectSolve(usingDirectSolve)
{
  mActiveMax = nobs/4 > 2 * hdim ? nobs/4 : 2 * hdim;

  mRhs        = SimpleVectorHandle(new SimpleVector(hdim));
  mYdReduced  = SimpleVectorHandle(new SimpleVector(hdim));
  mTraces     = SimpleVectorHandle(new SimpleVector(nobs));
  mYdReduced->setToZero();

  mIsSv       = new char[nobs];

  mP             = new double *[hdim];
  mP[0]          = new double[hdim*(hdim + 1)/2];

  mPfactored     = new double *[hdim];
  mPfactored[0]  = new double[hdim*(hdim + 1)/2];

  for( int j = 0; j < hdim; j++ ) {
    mP[j]          = &mP[0][j*(j+1)/2];
    mPfactored[j]  = &mPfactored[0][j*(j+1)/2];
  }
}


//////////////////////////////////////////////////////////////////////
SvmPrecond::~SvmPrecond()
{
  if( mP ) {
    delete [] mP[0];
    delete [] mP;
  }
  if( mPfactored ) {
    delete [] mPfactored[0];
    delete [] mPfactored;
  }
  delete [] mIsSv;
}


//////////////////////////////////////////////////////////////////////
void SvmPrecond::reset()
{
  mYdReduced->setToZero();

  for( int k = 0; k < (hdim*(hdim + 1))/2; k++ ) {
    mP[0][k] = 0.0; mPfactored[0][k] = 0.0;
  }
  for( int j = 0; j < hdim; j++ ) {
    mP[j][j] = 2.0;
  }
  for( int i = 0; i < nobs; i++ ) mIsSv[i] = kIsSvUnknown;

  mNsv           = 0;
  mGammaReduced  = 0.0;
  precon_count   = 0;
}


////////////////////////////////////////////////////////////////
void SvmPrecond::init(SvmMatrix & Y, SimpleVector & dinv, double mu)
{
  //double       StartUserTime, StartSystemTime;
  //double       EndUserTime, EndSystemTime;
  SvmMatrix::scalarT * M = Y.M();
  int * krowM = Y.krowM();
  int * jcolM = Y.jcolM();

  this->reset();
  //GetTime( &StartUserTime, &StartSystemTime);
  {
    double       **P      = this->getP();
    {
      SimpleVector & traces = *mTraces;
      for( int i = 0; i < nobs; i++ ) {
        double tracesi = 0;
        int kk = krowM[i];
        SvmMatrix::scalarT * row = M + kk;
        int * jcol = jcolM + kk;
        int nnz = krowM[i + 1] - krowM[i];

        for (int k = 0;  k < nnz;  k++) {
          int j = jcol[k];
          double temp = row[k] * row[k] * dinv[i];
          tracesi  += temp;
          P[j][j]  += temp;
        }
        traces[i]  = tracesi;
        traces[i] /= hdim;  // Scale down.
      }
    }
    this->chooseSupportVectors( *mTraces, mu );

    double ** Pcopy = this->getPworkSpace();
    const int binsize = 512;
    for( int i = 0; i < nobs; i++ ) {
      if( this->isSv(i) == kIsSv ) {
        int k = krowM[i];
        int nnz = krowM[i + 1] - k;
        spspr_1(hdim, dinv[i], M + k, jcolM + k, nnz, Pcopy[0]);
      }
      if( (i+1) % binsize == 0 || (i+1) == nobs ) {
        this->accumulateP();
      }
    }
  }  // end scope of dinv

  this->formReducedmYd( Y, dinv);

  //GetTime( &EndUserTime, &EndSystemTime );
  // if( gOoqpPrintLevel >= 100 ) {
  //   printf("Time to form preconditioner = %1.10g\n",
  //          EndUserTime - StartUserTime);
  // }

  this->factorP();
}


////////////////////////////////////////////////////////////////
void SvmPrecond::update(SvmMatrix & Y, SimpleVector   & dinv, double mu)
{
  //double       StartUserTime, StartSystemTime;
  //double       EndUserTime, EndSystemTime;

  this->chooseSupportVectors( *mTraces, mu );
  {
    double       **P      = this->getP();
    //GetTime( &StartUserTime, &StartSystemTime);
    int * krowM = Y.krowM();
    int * jcolM = Y.jcolM();
    SvmMatrix::scalarT * M = Y.M();

    for( int i = 0; i < nobs; i++ ) {
      if( this->isSv(i)  == kIsNewSv ) {
        int k = krowM[i];
        int nnz = krowM[i + 1] - k;
        spspr_1(hdim, dinv[i], M + k, jcolM + k, nnz, P[0]);
      } //end if  thisitioner->isSv(i)  == kIsNewSv
    } // end for i
  } //
  this->formReducedmYd( Y, dinv);

  //GetTime( &EndUserTime, &EndSystemTime );
  // if( gOoqpPrintLevel >= 100 ) {
  //   printf("Time to form updated preconditioner = %1.10g\n",
  //          EndUserTime - StartUserTime);
  // }
  this->factorP();
}


////////////////////////////////////////////////////////////////
void SvmPrecond::resetGamma(double mu)
{
  double mySvmtol = mu < 1 ? pow(mu, 0.5) : 1;
  SimpleVector & traces = *mTraces;

  //TESTCODE
  std::vector<double> tmp(nobs);
  for( int i = 0; i < nobs; i++)
    tmp[i] = traces[i];

  int my_max = mActiveMax > 0 ? mActiveMax : 1;
  if( my_max > nobs ) {
    my_max = nobs;
  }
  std::nth_element( tmp.begin(), tmp.end()-my_max, tmp.end());
  double nthmin = tmp[nobs-my_max];

  if(mActiveMax == 0) {
    nthmin *= 1.000001;
  }

  mScalesvmtol = .99999999 * nthmin/mySvmtol;
  //ENDTESTCODE
}



//////////////////////////////////////////////////////////////////////
void
SvmPrecond::chooseSupportVectors( SimpleVector & dinv, double mu )
{
  double mySvmtol = mu < 1 ? pow(mu, 0.5) : 1;

  double scale;
  if(mUsingDirectSolve) {
    scale = 0;
  } else {
    scale =  mScalesvmtol * mySvmtol;
  }

  for( int i = 0; i < nobs; i++ ) {
    if( dinv[i] >= scale ) {
      switch( mIsSv[i] ) {
        case kIsNotSv: mIsSv[i] = kIsNewSv; mNsv++; break;
        case kIsSv: break;
        case kIsSvUnknown:  mNsv++; mIsSv[i] = kIsSv; break;
        case kIsNewSv: /* was new, now old */ mIsSv[i] = kIsSv; break;
        default: assert( 0 && "Can't get here." ); break;
      } // end switch
    } else {
      mIsSv[i] = kIsNotSv;
      // end else
    }
  } // end for
  if( gOoqpPrintLevel >= 100 ) {
    printf("********************************\n");
    printf("Number of support vector = %d, s.v. tol = %1.6g,\n"
           "s.v. absolute = %1.6g\n",
           mNsv, scale, scale );
  }
}


//////////////////////////////////////////////////////////////////////
void
SvmPrecond::formReducedmYd( SvmMatrix & Y, SimpleVector & dinv)
{
  mYdReduced->setToZero();
  mGammaReduced = 0;
  int * krowM = Y.krowM();
  int * jcolM = Y.jcolM();
  SvmMatrix::scalarT * M =  Y.M();
  for(int i = 0; i < nobs; i++){
    switch(this->isSv(i)) {
      case kIsSv: case kIsNewSv:
        {
          int k = krowM[i];
          int nnz = krowM[i + 1] - k;
          spaxpy(dinv[i], M + k, jcolM + k, nnz, mYdReduced->elements());

          mGammaReduced += dinv[i];
        }
        break;
      case kIsNotSv: break;
      default:
        cout << "i: " << i << "isSv[i] " << this->isSv(i) << endl;
        assert( 0 && "Can't get here" );
        break;
    }
  }
}


//////////////////////////////////////////////////////////////////////
void SvmPrecond::accumulateP()
{
  for( int j = 0; j < hdim; j++ ) {
    double * Pjk = mP[j], *lPj = &mP[j][j];
    double * Pfactoredjk;
    for( Pfactoredjk = &mPfactored[j][0]; Pjk <= lPj; Pjk++, Pfactoredjk++ ) {
      *Pjk += *Pfactoredjk;
      *Pfactoredjk = 0;
    }
  }
}


//////////////////////////////////////////////////////////////////////
void SvmPrecond::factorP()
{
  //double       StartUserTime, StartSystemTime;
  //double       EndUserTime, EndSystemTime;

  char fortranUplo = 'U'; int info;
  int n = hdim*(hdim + 1)/2, ione = 1;
  //GetTime( &StartUserTime, &StartSystemTime);

  memcpy( mPfactored[0], mP[0], n * sizeof(double) );
  if( mNsv > 0 ) {
    double alpha = -1/mGammaReduced;

    dspr_(&fortranUplo, &hdim, &alpha, mYdReduced->elements(),
          &ione, mPfactored[0] );
    dpptrf_( &fortranUplo, &hdim, &mPfactored[0][0], &info );
  } else {
    if( gOoqpPrintLevel >= 100 ) {
      cout << "Preconditioner is diagonal!!!!\n";
    }
    for( int i = 0; i < hdim; i++ ) {
      // Since there are no support vectors, P is a diagonal matrix.
      mPfactored[i][i] = sqrt( mPfactored[i][i] );
    }
  }
  //GetTime( &EndUserTime, &EndSystemTime );
  // if( gOoqpPrintLevel >= 100 ) {
  //   printf("Time to factor preconditioner = %1.10g\n",
  //          EndUserTime - StartUserTime);
  // }
}


//////////////////////////////////////////////////////////////////////
void SvmPrecond::apply(SimpleVector &  xptr, SimpleVector & yptr)
{
  //double StartUserTime, StartSystemTime;
  //double EndUserTime, EndSystemTime;

  SimpleVector & rhs = *mRhs;

  //GetTime( &StartUserTime, &StartSystemTime);
  rhs.copyFrom( xptr );
  if( mNsv > 0 ) {
    char fortranUplo = 'U'; int info, one = 1;

    dpptrs_( &fortranUplo, &hdim, &one, mPfactored[0],
             &rhs[0], &hdim, &info);
  } else {
    for( int i = 0; i < hdim; i++ ) {
      rhs[i] /= (mPfactored[i][i] * mPfactored[i][i]);
    }
  }
  yptr.copyFrom( rhs );
  //GetTime( &EndUserTime, &EndSystemTime);
  // if( gOoqpPrintLevel >= 100 ) {
  //   printf("% 3d: Time to solve with preconditioner = %1.10g\n",
  //          precon_count, EndUserTime - StartUserTime);
  // }

  precon_count++;
}
extern int gOoqpPrintLevel;

//extern "C" void GetTime( double  * utime, double * stime );


////////////////////////////////////////////////////////////////
SvmLinearSolver::SvmLinearSolver(int m, int n, int usingDirectSolve):
  hyperplanedim(n), nobservations(m),
  mSvmmaxits(40),  mSolvestep(-1),
  mY(0), mYd(0), mDinv(0), mSvmPrecond(0),
  mScaleSvmTolMax(10000)
{
  mSvmPrecond
    = SvmPrecondHandle( new SvmPrecond( nobservations, hyperplanedim,
                                        usingDirectSolve ) );
  mX = SimpleVectorHandle( new SimpleVector(hyperplanedim) );
}


////////////////////////////////////////////////////////////////
SvmLinearSolver::~SvmLinearSolver()
{
}


////////////////////////////////////////////////////////////////
int SvmLinearSolver::calcMaxIts()
{
  return hyperplanedim/8 > 20 ? hyperplanedim/8 : 20;
}


////////////////////////////////////////////////////////////////
void SvmLinearSolver::newData(SvmMatrix * Y, SimpleVector * Yd,
                             SimpleVector * Dinv,
                             double gamma, double mu)
{
  // Display the number of outer iterations
  {
    static int mNumouterits = -1;
    if( gOoqpPrintLevel >= 100 )
      printf("My outer its = %d \n", ++mNumouterits);
  }
  SpReferTo(mY, Y);
  SpReferTo(mYd, Yd);
  SpReferTo(mDinv, Dinv);
  mGamma  = gamma;  mMu     = mu;

  mSolvestep = 0;
  //zero for predictor step, > 0 for corrector.

  if (gOoqpPrintLevel >= 100) {

    cout << "Scale svm tol " << mSvmPrecond->scaleSvmTol() << endl;
  }
  mSvmPrecond->init(*mY, *mDinv, mMu);
}


////////////////////////////////////////////////////////////////
void SvmLinearSolver::solve(SimpleVector& b)
{
  assert(mSvmPrecond);
  int ione = 1;

  const double rtolFactor = 1e-1;
  const double rtolMax    = 1e-1;
  const int    maxtries   = 5;
  const double atol       = 1e-12;
  double rtol = rtolMax < rtolFactor*mMu ? rtolMax : rtolFactor*mMu;
  double normb = dnrm2_(&hyperplanedim, b.elements(), &ione);
  SimpleVector & x = *mX;

  // If this is a corrector step, tighten the solution tolerance
  if( mSolvestep == 0 ) {
    // This is a predictor step; so reset the number of iterations
    mSvmmaxits = this->calcMaxIts();
  } else {
    // This is a corrector step so tighten the solution tolerance
    rtol *= 1e-2;
  }
  if ( gOoqpPrintLevel >= 100 ) {
    cout << ( mSolvestep > 0  ? "Corrector step\n" : "Predictor step\n");
    cout << "*****mKsprtol = " << rtol << endl;
  }
  for( int tries = 0; tries < maxtries; tries++ ) {
    if( mSolvestep > 0 || tries > 0 ) {
      // x contains a valid initial guess, use it
      //this->setInitialGuessNonzero();
    } else {
      x.setToZero();
    }
    if( gOoqpPrintLevel >= 100 ) {
      printf( "normb: %g rtol: %g atol: %g\n", normb, rtol, atol);
      printf( "Max its: %d\n", mSvmmaxits);
    }
    int cits = 0;
    int converged = 0;
    {
      double rprp = 1;
      SimpleVector r(hyperplanedim);
      SimpleVector p(hyperplanedim);
      SimpleVector pp(hyperplanedim);

      if (gOoqpPrintLevel >= 1000) {
        cout << "  KSP RHS " << normb << endl;
      }
      this->matMult(x, r);
      r.axpy(-1, b);
      while(1) {
        double normr = dnrm2_(&hyperplanedim, r.elements(), &ione);
        if (gOoqpPrintLevel >= 1000) {
          cout << "  KSP iteration " << cits << " normr " << normr << endl;
        }
        if (normr <= atol || normr <= rtol * normb) {
          converged = 1;
          break;
        }
        cits++;
        if (cits > mSvmmaxits)
          break;
        mSvmPrecond->apply(r, p);
        double rp = r.dotProductWith(p);
        if (cits > 1) {  /* there are prior directions, i.e. pp in nonzero. */
          double beta = rp/rprp;
          p.axpy(beta, pp);
        }
        this->matMult(p, pp);
        double pAp = p.dotProductWith(pp);
        double alpha = -rp/pAp;
        x.axpy(alpha, p);
        r.axpy(alpha, pp);
        /* Set up for next iteration */
        rprp = rp;
        pp.copyFrom(p);
      }
    }
    mSvmmaxits -= cits;
    if( mSvmmaxits < 1 ) mSvmmaxits = 1;

    if( converged ) {
      break;
    } else if( mSvmPrecond->nSv() < nobservations ) {
      int nsv =  mSvmPrecond->nSv();
      int newActiveMax = nsv +  hyperplanedim/2;
      if( newActiveMax <= nsv ) {
        newActiveMax = nsv + 1;
      }
      mSvmPrecond->setActiveMax(newActiveMax);
      mSvmPrecond->resetGamma(mMu);
      mSvmPrecond->update(*mY, *mDinv, mMu);

      mSvmmaxits = this->calcMaxIts();
    } else { // no hope
      fprintf(stderr, "Failed to solve preconditioning with matrix itself."
              "  Aborting ...\n");
      exit( 1 );
    } // end else no hope
  } // end for tries < maxtries
  b.copyFrom(x);
  mSolvestep++;
}


////////////////////////////////////////////////////////////////
int SvmLinearSolver::converged()
{
  int reason = 1;

  return reason > 0;
}


////////////////////////////////////////////////////////////////
void SvmLinearSolver::matMult(SimpleVector & avec, SimpleVector & yvec)
{
  const int binsize = 8 * 1024;

  /* remember mY is stored in transposed format */
  /* y =  mY'*W*mY*a - (1/gamma)*mYd*mYd'*a + 2a */
  int nobs = nobservations, hdim = hyperplanedim;

  SimpleVector & Yd   = *mYd;
  SimpleVector & dinv = *mDinv;

  int * krowM = mY->krowM();
  int * jcolM = mY->jcolM();
  SvmMatrix::scalarT * M = mY->M();

  {
    SimpleVectorHandle hyvec_copy( new SimpleVector( hdim ) );
    SimpleVector & yvec_copy = * hyvec_copy;

    yvec.copyFrom( avec );
    yvec.scale(2.0);

    yvec_copy.setToZero();
    int j;
    SimpleVector Yj(hdim);
    for (j = 0;  j < nobs;  j++) {
      int k = krowM[j];
      SvmMatrix::scalarT * row = M + k;
      int * jcol = jcolM + k;
      int nnz = krowM[j + 1] - krowM[j];
      double prod = spdot(avec.elements(), row, jcol, nnz);
      prod *= dinv[j];
      spaxpy(prod, row, jcol, nnz, yvec_copy.elements());
      if (j + 1 % binsize == 0) {
        yvec.axpy(1.0, yvec_copy);
        yvec_copy.setToZero();
      }
    }
    if (j % binsize != 0.0) {
      yvec.axpy(1.0, yvec_copy);
    }
  } // end scopy of yvec_copy

  {
    double Ydta = Yd.dotProductWith( avec );
    Ydta = -Ydta/mGamma;
    yvec.axpy( Ydta, Yd );
  } // end scope of Ydta
}
