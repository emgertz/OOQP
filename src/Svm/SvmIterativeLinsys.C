/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include <cmath>
#include <iostream>

using namespace std;

#include "SvmIterativeLinsys.h"
#include "DenseSymMatrix.h"
#include "DenseGenMatrix.h"
#include "SvmData.h"
#include "SvmVars.h"
#include "SvmResiduals.h"
#include "SimpleVector.h"
#include "DoubleMatrix.h"
#include "SvmLinearSolver.h"

extern "C" {
  void dsyr_( char *, int *, double *, double *, int *, double *, int *);
}


SvmIterativeLinsys::SvmIterativeLinsys(SvmData *prob, int usingDirectSolve) :
  mYd(0), mGamma(0), mDinv(0)
{
  int hyperplanedim = prob->hyperplanedim;
  int nobservations = prob->nobservations;

  // allocates the messy diagonal matrix that arises in the block
  // elimination
  mDinv = SimpleVectorHandle( new SimpleVector( nobservations ) );
  mYd   = SimpleVectorHandle( new SimpleVector( hyperplanedim ) );

  mSolver = SvmLinearSolverHandle( new SvmLinearSolver( nobservations,
							hyperplanedim, usingDirectSolve) );
}


SvmIterativeLinsys::~SvmIterativeLinsys()
{
}


void SvmIterativeLinsys::factor(Data *prob_in, Variables *vars_in)
{
  // does a Cholesky factorization of the matrix (X^T D^{-1} X)
  char fortranUplo = 'U'; int info;
  SvmData *prob = (SvmData *) prob_in;
  SvmVars *vars = (SvmVars *) vars_in;
  int i, j;

  // compute and store mDiv, using vars

  double *ps = vars->s->elements(), *pv = vars->v->elements();
  double *pz = vars->z->elements(), *pu = vars->u->elements();
  {
    double *pDinv = mDinv->elements(), *lDinv = pDinv + mDinv->length();
    for( ; pDinv < lDinv; pDinv++, ps++, pv++, pz++, pu++ ) { 
      *pDinv = 1 / ( *ps / *pv + *pz / *pu );
    }
  }
  // now compute lower triangle of X^T D^{-1} X, storing in L
  prob->XTransMult(0.0, *mYd, 1.0, *mDinv);
  SimpleVector & dinv = *mDinv;
  mGamma = 0.0;
  for( int i = 0; i < dinv.length(); i++ ) {
    mGamma += dinv[i];
  }
  mSolver->newData( prob->mY, mYd, mDinv, mGamma, vars->mu());
}

void SvmIterativeLinsys::solve(Data *prob_in, Variables *vars_in, 
		      Residuals *resids_in, Variables *step_in)
{
  SvmData *prob = (SvmData *) prob_in;
  SvmResiduals *resids = (SvmResiduals *) resids_in;
  SvmVars *vars = (SvmVars *) vars_in;
  SvmVars *step = (SvmVars *) step_in;

  // So we don't go mad typing, make local references to all interesting
  // quantities.
  SimpleVector & Dinv   = *mDinv;

  SimpleVector & uVars  = *vars->u;
  SimpleVector & vVars  = *vars->v;
  SimpleVector & sVars  = *vars->s;
  SimpleVector & zVars  = *vars->z;

  SimpleVector & wStep  = *step->w;
  SimpleVector & uStep  = *step->u;
  SimpleVector & vStep  = *step->v;
  SimpleVector & sStep  = *step->s;
  SimpleVector & zStep  = *step->z;
  double   &       betaStep = step->beta;

  SimpleVector & wRes  = *resids->wRes;
  SimpleVector & rUZ   = *resids->rUZ;
  SimpleVector & rSV   = *resids->rSV;
  SimpleVector & sRes  = *resids->sRes;
  SimpleVector & zRes  = *resids->zRes;
  double           betaRes = resids->betaRes;

  { // context UinvZ
    SimpleVector & UinvZ  = sStep;
    // UinvZ = U\inv * Z
    UinvZ.copyFrom( zVars );
    UinvZ.componentDiv( uVars );    
    { // context for rzhat
      SimpleVector & rzhat  = zStep;
      // rzhat = zRes + Z\inv * rUZ
      rzhat.copyFrom( zRes );
      rzhat.axdzpy( 1.0, rUZ, zVars );

      { // context for omegaRes
	SimpleVector & omegaRes = uStep;
	// romega = D\inv ( sRes + V\inv * vRes + U\inv * Z * rzhat )
	omegaRes.copyFrom( sRes );
	omegaRes.axdzpy( 1.0, rSV, vVars );
	omegaRes.axzpy ( -1.0, UinvZ, rzhat );
	omegaRes.componentMult( Dinv );
	
	wStep.copyFrom( wRes );
	prob->YTransMult( 1.0, wStep, 1.0, omegaRes );
	
	double rhatbeta;
	rhatbeta = betaRes - prob->dotCategories( omegaRes );
	wStep.axpy( rhatbeta/mGamma, *mYd );
	
	mSolver->solve(wStep);
	
	// betaStep = (rhatbeta + < wStep, Yd > ) / mGamma
	betaStep = (rhatbeta + wStep.dotProductWith( *mYd )) / mGamma;

	// Finish the computation of vStep
	// vStep = Omega\inv( categories * betaStep - Y * wStep)
	//              + omegaRes
	vStep.copyFrom( *prob->categories );
	vStep.scale( betaStep );
	prob->YMult( 1.0, vStep, -1.0, wStep );
	vStep.componentMult( Dinv );
	vStep.axpy( 1.0, omegaRes );
      } // omegaRes dissappears, freeing uStep
    } // rzhat disappears, and its value is in zStep
    // Finish the computation of zStep
    // zStep = - U\inv Z ( rzhat + vStep )
    zStep.axpy( 1.0, vStep );
    zStep.componentMult( UinvZ );
  } // UinvZ dissappears (releasing sStep)

  // sStep = (rhs.rSV - vStep * sVars) / vVars
  sStep.copyFrom( rSV );
  sStep.axzpy( -1.0, vStep, sVars );
  sStep.componentDiv( vVars );
  
  // uStep = (rhs.rUZ - zStep * uVars) / zVars
  uStep.copyFrom( rUZ );
  uStep.axzpy( -1.0, zStep, uVars );
  uStep.componentDiv( zVars );
}






