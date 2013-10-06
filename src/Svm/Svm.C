/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "Svm.h"
#include "SvmData.h"
#include "SvmResiduals.h"
#include "SvmVars.h"
#include "SvmLinsys.h"
#include "SvmIterativeLinsys.h"
#include "DoubleMatrix.h"


// SVMFactory::SVMFactory( )
// {
// }

Data * Svm::makeData(int hyperplanedim, int nobservations, int nnz, double penalty )
{
  return new SvmData(hyperplanedim, nobservations, nnz, penalty );
}

// Data * Svm::makeData(int hyperplanedim, int nobservations,
// 		     double * X, double * d, double penalty )
// {
//   return new SvmData(hyperplanedim, nobservations, X, d, penalty );
// }

Data * Svm::makeRandomData(int hyperplanedim, int nobservations,
			   double penalty )
{
  SvmData * data = (SvmData *) this->makeData(hyperplanedim, nobservations, 
					      hyperplanedim * nobservations,
					      penalty);
  data->datarandom();
  return data;
}

Data * Svm::makeDataFromText(char  filename[], double penalty, int& iErr,
                             bool dense_input)
{
  if (dense_input) {
    return SvmData::denseTextInput(filename, penalty, iErr);
  } else {
    return SvmData::textInput(filename, penalty, iErr);
  }
} 

Residuals * Svm::makeResiduals( Data * prob_in )
{
  SvmData * prob = (SvmData *) prob_in;

  return new SvmResiduals( prob->hyperplanedim, prob->nobservations );
}

Variables * Svm::makeVariables( Data * prob_in )
{
  SvmData * prob = (SvmData *) prob_in;
  return new SvmVars( prob->hyperplanedim, prob->nobservations  );
}

Variables * Svm::makeVariables( Data * prob_in, 
				double w[], double v[], 
				double z[], double u[],
				double s[] )
{
  SvmData * prob = (SvmData *) prob_in;
  return new SvmVars( prob->hyperplanedim, prob->nobservations,
		      w, v, z, u, s  );
}


Svm::~Svm()
{
}

LinearSystem * SvmDirect::makeLinsys( Data * prob_in )
{
  SvmData & prob = dynamic_cast<SvmData&>(*prob_in);
  
  return new SvmLinsys(&prob);
}


LinearSystem * SvmIterative::makeLinsys( Data * prob_in )
{
  SvmData & prob = dynamic_cast<SvmData&>(*prob_in);
  
  return new SvmIterativeLinsys(&prob, mUsingDirectSolve);
}


void SvmStartStrategy::doIt( Solver * /* solver */,
			     ProblemFormulation * /* formulation */, 
			     Variables * iterate, Data * /* prob*/,
			     Residuals * /* resid */, Variables * /*step */)
{
  iterate->interiorPoint(2.0, 2.0);
};
