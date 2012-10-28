/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef SVMITERATIVELINSYS_H
#define SVMITERATIVELINSYS_H

#include "LinearSystem.h"
#include "SimpleVectorHandle.h"
#include "SvmLinearSolverHandle.h"

class Data;
class Variables;
class Residuals;
class DoubleLinearSolver;

/**
 * @ingroup Svm
 *
 * LinearSystem class for Svm.
 *
 */

#include "DenseSymMatrixHandle.h"

class SvmData;
  
class SvmIterativeLinsys : public LinearSystem
{
public:
  SvmIterativeLinsys(SvmData * prob, int usingDirectSolve);
  ~SvmIterativeLinsys();
  void factor(Data *prob, Variables *vars);

  void solve(Data *prob, Variables *vars, Residuals *rhs,
	     Variables *step);

private:
  /** */
  SimpleVectorHandle mYd;
  
  /** */
  double mGamma;

  /** pointer to dense symmetric positive definite solver used to
      factor the compressed matrix */
  SvmLinearSolverHandle mSolver;

  /** stores the complicated diagonal matrix that arises during the
   *  block elimination */
  SimpleVectorHandle mDinv;
};
  

#endif
