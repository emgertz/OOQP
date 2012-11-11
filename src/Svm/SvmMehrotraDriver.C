/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include <cstring>
#include <iostream>
#include <fstream>
#include <cstdio> 
#include <cstdlib>
using namespace std;

#include "SvmData.h"
#include "SvmVars.h"
#include "SvmResiduals.h"
#include "SvmLinsys.h"
#include "MehrotraSolver.h"
#include "Svm.h"
#include "SimpleVector.h"

//extern "C" void GetTime( double  * utime, double * stime );

extern int DenseStorageInstances; 
extern int gOoqpPrintLevel;

int main( int argc, char *argv[] )
{
  Svm     * svm         = 0;
  SvmData * prob        = 0;
  int quiet = 0, print_soln = 0, using_direct_solve = 0, dense_input = 0;
  
  double StartUserTime, StartSystemTime;
  double EndUserTime, EndSystemTime;
  
  char    * outfilename = 0;
  int argsOk = 1;
  {
    int iarg = 1;
    
    while (iarg < argc && argv[iarg][0] == '-') {
      // it is a option. Check against recognized options
      if( 0 == strcmp( argv[iarg], "-quiet" ) ||
	  0 == strcmp( argv[iarg], "--quiet" ) ) {
	quiet = 1;
      } else if ( 0 == strcmp( argv[iarg], "-print-level" ) ||
		  0 == strcmp( argv[iarg], "--print-level" ) ) {
	iarg++;
	gOoqpPrintLevel = atoi(argv[iarg] );
      } else if ( 0 == strcmp( argv[iarg], "-print-solution" ) ||
		  0 == strcmp( argv[iarg], "--print-solution" ) ) {
        
	print_soln = 1;
        
      } else if ( 0 == strcmp( argv[iarg], "-direct-solve" ) ||
		  0 == strcmp( argv[iarg], "--direct-solve" ) ) {
	using_direct_solve = 1;
      } else if ( 0 == strcmp( argv[iarg], "-dense-input" ) ||
                  0 == strcmp( argv[iarg], "--dense-input" ) ) {
        dense_input = 1;
      } else {
	cerr << argv[0] << ": "
	     << argv[iarg] << " is not a recognized option.\n";
	argsOk = 0;
      }
      iarg++;
    }
    if( iarg >= argc ) argsOk = 0;  // Not enough arguments 
    if( argsOk ) { // the options were read successfully
      if( 0 == strcmp( argv[iarg] , "random" ) ) {
	int hyperplanedim, nobservations;

	if( iarg + 2 >= argc ) { 
	  // Not enough arguments
	  argsOk = 0;
	} else { // there are enough args
	  hyperplanedim = atoi(argv[iarg + 1]); 
	  nobservations = atoi(argv[iarg + 2]);
	  if(hyperplanedim <= 0) {
	    cerr << " Hyperplane dimension must be positive;"
		 << " hyperplanedim=" <<  hyperplanedim << endl;
	    argsOk = 0;
	  }
	  if(nobservations <= 1) {
	    cerr << " Must have at least two observations;"
		 << " nobservations=" <<  nobservations << endl;
	    argsOk = 0;
	  }
	  if(nobservations < hyperplanedim) {
	    cerr << " Number of observations must be no smaller"
		 << " than hyperplane dimension\n";
	    argsOk = 0;
	  }
	  if( argsOk ) {
	    prob = (SvmData *) svm->makeRandomData(hyperplanedim,
						   nobservations, 1.0);
	  }
	} // there are enough args
      } else { // The data is to be read from a file
	// There can be at most one more argument
	if( argc > iarg + 2 ) { // Too many args
	  argsOk = 0;
	} else { // We have the right number of arguments
	  double penalty  = 1.0;
	  if( iarg + 1 < argc ) { // the penalty parameter was specified
	    char * endptr;
	    char * penstr = argv[iarg+1];
	    penalty = strtod( penstr, &endptr );
	    if( penstr == endptr || // No conversion could be done or
		*endptr != '\0' ) { // We didn't use the whole argument
	      cerr << "Error: I couldn't parse " << penstr 
		   << "as a floating point number\n";
	      argsOk = 0;
	    }
	  } // end if the penalty parameter was specified
	  if( argsOk ) { // syntax of the args is good
            svm =  new SvmIterative(using_direct_solve);
	    char * filename = argv[iarg];
	    // Get a name for the output file.
	    outfilename = new char[strlen(filename) + 5]; 
	    strcpy( outfilename, filename );
	    strcat( outfilename, ".out" );
	    // Try to read the input file
	    int iErr;
	    // Note that we, probably incorrectly, use twice the penalty
	    // parameter internally.
	    prob = (SvmData *) 
	      svm->makeDataFromText(filename, 2.0 * penalty, iErr, dense_input);
	    if(iErr != svminputok) {
	      cerr << " Error reading input file " << filename 
		   <<": TERMINATE\n";
	      return 1; // The args parse, but we can't read the file
	    }	  
	  } // if syntax of the args is good
	} // end else have the right number of arguments
      }  // end else the data is to be read from a file
    } // end if the options were read successfully
  } // end of the scope of iarg
  if( !argsOk ) {
    cerr << "\nUsage: \n\n";
    cerr << "    " << argv[0] << " [ --quiet ] [ --print-solution ]"
	 << " filename [ penalty ] \n\nor\n\n";
    cerr << "    " << argv[0] << " [ --quiet ] [ --print-solution ] "
	 << "random hdim nobs\n\n";
    cerr << "where \"random\" is a literal keyword.\n\n";

    delete svm;
    delete prob;

    return 1;
  }

  //GetTime( &StartUserTime, &StartSystemTime );

  MehrotraSolver * s     = new MehrotraSolver( svm, prob );
  SvmStartStrategy startStrategy;
  //s->useStartStrategy(&startStrategy);

  SvmVars       * vars  = (SvmVars *) svm->makeVariables( prob );
  Residuals     * resid = svm->makeResiduals( prob );
  if( !quiet ) {
    s->monitorSelf();
  }
  s->setMuTol(1e-6);
  s->setArTol(1e-6);
  int status = s->solve(prob, vars, resid);
  
  //GetTime( &EndUserTime, &EndSystemTime );

  // print the interesting variables
  if( (!quiet && vars->hyperplanedim < 20) ||
      !outfilename || print_soln ) {
    cout.precision(4);
    vars->printCoefs();
  }
  if( outfilename ) {
    {
      ofstream outfile( outfilename );
      outfile.precision(16);
      outfile << vars->hyperplanedim << endl;
      vars->w->writeToStream( outfile );
      outfile << vars->beta << endl;
    }
    delete [] outfilename;
  } 

  delete vars;  
  delete resid;
  delete s;
  delete prob;
  delete svm;

  // printf("Time to solve QP = %1.8g user %1.8g system\n",
  //        EndUserTime - StartUserTime, EndSystemTime - StartSystemTime);

  exit(status);
  return status;
}


