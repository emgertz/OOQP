#ifndef SVMPRECONDITIONER
#define SVMPRECONDITIONER

#include "SimpleVectorHandle.h"
#include "SvmLinearSolverHandle.h"
#include "SvmMatrix.h"

enum { kIsNotSv = 0, kIsSv = 1, kIsNewSv = 2, kIsSvUnknown = -1 };


class SvmPrecond : public IotrRefCount {
public:
  void reset();
  void accumulateP();
  void factorP();
  void formReducedmYd( SvmMatrix & pY, SimpleVector & dinv);
  virtual void init(SvmMatrix & Y, SimpleVector & dinv, double mu);
  virtual void update(SvmMatrix & Y, SimpleVector & dinv, double mu);

  void apply( SimpleVector & x, SimpleVector & y );
  double ** getP() { return mP; }
  double ** getPworkSpace() { return mPfactored; }
  double    gammaReduced() { return mGammaReduced; }
  int isSv( int index ) { return (int) mIsSv[index]; }
  int nSv() { return mNsv; }
  void resetGamma(double mu);
  void chooseSupportVectors(SimpleVector & dinv, double mu);
  int activeMax() { return mActiveMax; }
  void setActiveMax( int activeMax ) { 
    mActiveMax = activeMax > nobs ? nobs : activeMax;
  }
  double scaleSvmTol() { return mScalesvmtol; }
  void setScaleSvmTol( double scale ) { mScalesvmtol = scale; }
  SvmPrecond( int m, int n, int usingDirectSolve );
  ~SvmPrecond();
protected:
  double **mP;
  double **mPfactored;
  SimpleVectorHandle mYdReduced;
  SimpleVectorHandle mRhs;
  SimpleVectorHandle    mTraces;
  double mGammaReduced;
  int    nobs;
  int    hdim;
  char   * mIsSv;
  static double  mScalesvmtol;
  int mActiveMax;
  int mNsv;
  int mUsingDirectSolve;
};
typedef SmartPointer<SvmPrecond> SvmPrecondHandle;

/**
 * A linear solver that uses the preconditioned conjugate gradient method in a
 * fashion that is specific to solving the Support Vector Machine optimization
 * problem */
class SvmLinearSolver : public IotrRefCount
{
 public:
  SvmLinearSolver(int m, int n, int usingDirectSolve);
  virtual int newData(SvmMatrix * Y,
		      SimpleVector *  Yd, SimpleVector *
		      Dinv, double gamma, double mu);
  virtual void solve(SimpleVector& b);
  virtual int converged();
  virtual int calcMaxIts();

  void matMult(SimpleVector & a, SimpleVector & y);

  virtual ~SvmLinearSolver();

 protected:
  int hyperplanedim, nobservations;
  int mSvmmaxits;
  int mSolvestep;  // mSolvestep = 0 if on predictor, = 1 corrector.
  /** The number of support vectors */
  double mScaleSvmTolMax;
  /** The arrays that define the matrix */
  SvmMatrixHandle  mY;
  SimpleVectorHandle    mYd, mDinv;
  double mGamma;
  double mMu;
  SvmPrecondHandle mSvmPrecond;
  SimpleVectorHandle mX;
};

#endif

