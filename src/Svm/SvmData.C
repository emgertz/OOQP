/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "SvmData.h"
// #include "math.h"
#include <cstring>
#include <cstdio>
#include <cerrno>
#include <iostream>
#include <fstream>
using namespace std;

#include <cassert>
#include <stdlib.h>
// include for IO purposes
#include <fstream>

#include "SparseGenMatrixT.h"
#include "SimpleVector.h"
#include "LinearAlgebraPackage.h"   

//define a constructor

SvmData::SvmData(int hyperplanedim_in, int nobservations_in, int nnz,
		 double penalty /* = 1.0 */ )
{
  hyperplanedim = hyperplanedim_in;
  nobservations = nobservations_in;
  mPenalty      = penalty;

  assert ( nobservations_in >= hyperplanedim_in);
  assert ( nobservations_in >= 2);
  assert ( hyperplanedim_in >= 1);

  mY = SvmMatrixHandle( new SvmMatrix(nobservations, hyperplanedim, nnz) );

  categories = SimpleVectorHandle( new SimpleVector( nobservations ) );
}

/** s_fgets -- a reimplementation of fgets that checks for long lines 
 *
 *  size - the size of buffer s.  
 *
 * At most size - 1 characters will be
 *  read, and the buffer will always be terminated by '\O'.
 */
static 
char * s_fgets(char * s, int size, FILE * stream, int& err)
{
  err = 0;
  int c = '\0';
  int i = 0;
  while (i < size - 1 && EOF != (c = getc(stream))) {
    s[i++] = c;
    if (c == '\n')
      break;
  }
  s[i] = '\0';
  if (i == 0)
    return NULL;
  else if (c == '\n' || c == EOF)
    return s;
  else {
    err = -1;
    return NULL;
  }
}

static int parseInt(char * str, int & status)
{
  char * endptr;
  errno = 0;
  int result = strtol(str, &endptr, 10);
  if (errno == 0 && *endptr == 0)
    status = 0;
  else
    status = 1;
  return result;
}

static double parseDouble(char * str, int & status)
{
  char * endptr;
  errno = 0;
  int result = strtol(str, &endptr, 10);
  if (errno == 0 && *endptr == 0)
    status = 0;
  else
    status = 1;
  return result;
}

/* The textinput routine is called when we read input from a file
 * 
 * The file is/should be the same (sparse) format as SVMlight.  The
 * largest line permitted is 16kb. */

SvmData * SvmData::textInput( char filename[], double penalty, int& iErr)
{
  int i;
  FILE * file;
  SvmData * result;
  char buffer[16 * 1024];
  char * linetok;
  char * itemtok;
  int status;

  // assume that the input will be OK
  iErr = svminputok;

  //open the input file; bomb out if it can't be found
  file = fopen( filename, "r");
  if(!file) {
    fprintf( stderr, " textInput: Error reading %s: fopen failed\n", filename);
    iErr = svmfileopenerror;
    return 0;
  }
  // scan to find the sparsity and size
  int nnz = 0;
  int nobservations = 0;
  int hyperplanedim = 0;
  int lineno = 0;
  while( s_fgets(buffer, sizeof(buffer), file, status) ) {
    ++lineno;
    nobservations ++;
    char * tok;
    // ignore the first, it is the category
    tok = strtok_r(buffer, " \t\n", &linetok);
    if (!tok) {
      fprintf( stderr, " textInput: Error reading %s on line %d\n", 
	       filename, lineno );
      iErr = svmfileinputerror;
      return 0;
    }

    while ((tok = strtok_r(NULL, " \t\n", &linetok))) {
      nnz ++;
      status = 1;
      char * cindex = strtok_r(tok, ":", &itemtok);
      
      if (cindex) {
	int index = parseInt(cindex, status);
	if (index > hyperplanedim) {
	  hyperplanedim = index;
	}
      }
      if (status != 0) {
	fprintf( stderr, " textInput: Error reading %s on line %d\n", 
		 filename, lineno );
	iErr = svmfileinputerror;
	return 0;
      }
    }
  }
  if (status != 0) {
    fprintf( stderr, " textInput: Error reading %s on line %d\n", 
	     filename, lineno );
    iErr = svmfileinputerror;
    return 0;
  }
  result = new SvmData(hyperplanedim, nobservations, nnz, penalty );          
  SimpleVector & categories = *result->categories;
  SimpleVectorHandle hTempRowY( new SimpleVector( result->hyperplanedim ) );
  SimpleVector & tempRowY = *hTempRowY;

  if (0 != fseek(file, 0L, SEEK_SET)) {
	fprintf( stderr, " textInput: Error rewinding %s\n", filename );
	iErr = svmfileinputerror;
	return 0;
  }    
  i = 0;
  lineno = 0;
  while(s_fgets(buffer, sizeof(buffer), file, status)) {
    ++lineno;
    tempRowY.setToZero();
    status = 1;
    char * tok = strtok_r(buffer, " \n\t", &linetok);
    if (tok) {
      int cat = parseInt(tok, status);
      if (status == 0) {
	if (cat == 1 || cat == -1) {
	  categories[i] = cat;
	} else {
	  fprintf( stderr, " textInput: Label not -1 or 1 in %s on line %d\n", 
		   filename, lineno );
	  iErr = svmlabelerror;
	  return 0;
	}
      }
    }
    if (0 != status) {
      fprintf( stderr, " textInput: Error reading %s on line %d\n", 
	       filename, lineno );
      iErr = svmfileinputerror;
      return 0;
    }
    while((tok = strtok_r(NULL, " \n\t", &linetok))) {
      status = 1;
      char * cindex = strtok_r(tok, ":", &itemtok);
      char * cvalue = strtok_r(NULL, ":", &itemtok);
      if (cindex && cvalue) {
	int index = parseInt(cindex, status);
	if (0 == status) 
	  tempRowY[index - 1] = parseDouble(cvalue, status);
      }
      if (0 != status) {
	fprintf( stderr, " textInput: Error reading %s on line %d\n", 
		 filename, lineno );
	iErr = svmfileinputerror;
	return 0;
      }
    }
    result->mY->atPutDense(i, 0, tempRowY.elements(), 1, 1, hyperplanedim);
    i++;
  }
  if (status != 0) {
    fprintf( stderr, " textInput: Error reading %s on line %d\n", 
	     filename, lineno );
    iErr = svmfileinputerror;
    return 0;
  }
  if (0 != fclose(file)) {
    fprintf( stderr, " textInput: Error closing %s\n", filename );
    iErr = svmfileinputerror;
    return 0;
  }

  return result;
}

/* The denseTextinput routine is called when we read input from a file
* 
* The file must start with two integers that denote
*
* - nobservations : number of observations
* - hyperplanedim : hyperplane dimension (that is, the size of each vector)
*
* must satisfy nobservations >= 2, hyperplanedim >= 1, nobservations
* >= hyperplanedim.
*
* Then for each observation one must list the "hyperplanedim" entries
* in each vector, followed by a label for that vector. The labels
* must take two distinct values (not necessarily +1 and -1).  */

SvmData * SvmData::denseTextInput( char filename[], double penalty, int& iErr)
{
    int i;
    int nobs_file, ndim_file;
    double label1, label2;
    FILE * file;
    SvmData * result;
    
    // assume that the input will be OK
    iErr = svminputok;
    
    //open the input file; bomb out if it can't be found
    file = fopen( filename, "r");
    if(!file) {
        fprintf( stderr, " textInput: Error reading %s: fopen failed\n", filename);
        iErr = svmfileopenerror;
        return 0;
    }
    
    // read dimensions : first the number of observations
    if ( fscanf(file, "%d", &nobs_file) != 1) {
        fprintf( stderr, " textInput: Error reading %s:"
                 " couldn't find nobservations\n",  filename);
        iErr = svmfileinputerror;
        return 0;
    }
    if ( fscanf(file, "%d", &ndim_file) != 1) {
        fprintf( stderr, " textInput: Error reading %s:"
                 " couldn't find hyperplanedim\n", filename);
        iErr = svmfileinputerror;
        return 0;
    }
    if(nobs_file <= 1) {
        fprintf( stderr, " Error reading %s:"
                 " Must have at least two observations\n", filename);
        fprintf( stderr, " nobservations=%d\n", nobs_file);
        iErr = svmfileinputerror;
        return 0;
    }
    if(ndim_file <= 0) {
        fprintf( stderr, " Error reading %s:"
                 " Data space dimension must be at least 1\n", filename);
        fprintf( stderr, " hyperplanedim=%d\n", ndim_file);
        iErr = svmfileinputerror;
        return 0;
    }
    // create a data structure with these dimensions
    result = new SvmData(ndim_file, nobs_file, nobs_file * ndim_file, penalty );          
    SimpleVector & categories = *result->categories;
    
    // and temporary storage for each row of the matrix;
    { 
        SimpleVectorHandle hTempRowY( new SimpleVector( result->hyperplanedim ) );
        SimpleVector & tempRowY = *hTempRowY;
        
        // read the matrix and labels
        for(i=0; i<result->nobservations; i++) {
            
            // read a row
            for(int j=0; j<result->hyperplanedim; j++) {
                if(fscanf(file, "%le", &tempRowY[j]) != 1) {
                    fprintf( stderr, " Error reading %s at observation %d\n",
                             filename, i);
                    fprintf( stderr, " Apparently too little data in the file\n");
                    iErr = svmfileinputerror;
                    return 0;
                }
            }
            // stuff the row into a column of Yt
            result->mY->atPutDense(i, 0, &tempRowY[0], 1, 1, result->hyperplanedim);
            
            //read the right-hand side element
            {
                double tlabel;
                if(fscanf(file, "%le", &tlabel) != 1) {
                    fprintf( stderr,
                             " Error reading %s at observation %d\n", filename, i);
                    fprintf( stderr, " Apparently too little data in the file\n");
                    iErr = svmfileinputerror;
                    return 0;
                }
                categories[i] = tlabel;
            }
        } // end for "read the matrix
  } // end scope of tempRowY
    fclose(file);
    
    // now check that there are just two distinct labels;
    
    label1 = categories[0];
    for (i = 0;
	 i < result->nobservations && categories[i] == label1;
	 i++ ) 
        ;
    
    if(i < result->nobservations) 
        label2 = categories[i];
    else {
        // only found one label
        fprintf(stderr, " Error reading %s: Found only one label: %g\n",
                filename, label1);
        iErr = svmlabelerror;
        return 0;
    }
    
    for(i=0; i<result->nobservations; i++) {
        if(categories[i] != label1 && categories[i] != label2) {
            // found a third label
            fprintf( stderr, " Error reading %s: "
                     "Found more than two labels: %g, %g, %g,...\n", 
                     filename, label1, label2, categories[i]);
            iErr = svmlabelerror;
            return 0;
        }
    }
    
    //printf(" Found exactly two labels: %g, %g\n", label1, label2);
    
    // Now, make the two labels be exactly 1, and -1
    // Make the greater of the two labels be 1 and the lesser -1.
    int negative_label = (label1 < label2) ? label1 : label2;

    for(i=0; i<result->nobservations; i++) {
        if(categories[i] == negative_label) {
            categories[i] = -1.0;
        } else {
            categories[i] =  1.0;
        }
    }
    return result;
}


SvmData::~SvmData()
{
}

// calculate the norm of the data for the SvmData class

double SvmData::datanorm()
{
  return mY->abmaxnorm();
}

  
void SvmData::YMult( double beta, SimpleVector& y,
		     double alpha, SimpleVector& x )
{
  SimpleVectorHandle svtemp( new SimpleVector(nobservations) );
  SimpleVector & temp = *svtemp;

  mY->mult( 0, temp, alpha, x );
  temp.componentMult( *categories );

  int i;
  for( i = 0; i < nobservations; i++ ) {
    y[i] = temp[i] + beta * y[i];
  }
}


void SvmData::XTransMult( double beta, SimpleVector  & y,
			  double alpha, SimpleVector & x )
{
  mY->transMult( beta, y, alpha, x );
}


void SvmData::YTransMult( double beta, SimpleVector& y,
			  double alpha, SimpleVector& x )
{
  // Probably resonable here to use a temp.
  SimpleVectorHandle svtemp( new SimpleVector( nobservations) );
  SimpleVector & temp = *svtemp;

  temp.copyFrom(x);
  temp.componentMult(*categories);

  mY->transMult( beta, y, alpha, temp );
}

void SvmData::datarandom()
{
  double drand(double *), ix;
  int i;

  ix = 89176823.0;

  // fill out the matrix Y with random numbers
  mY->randomize( -1.0, 1.0, &ix );
  categories->randomize( -1.0, 1.0, &ix );
  // set elements of the "categories" to -1 or 1
  double * pcat = &(*categories)[0];
  for(i=0; i<nobservations; i++) {
    if( pcat[i] < 0 ) {
      pcat[i] = -1.0;
    } else {
      pcat[i] =  1.0;
    }
  }
}

double SvmData::dotCategories( SimpleVector & v )
{
  return categories->dotProductWith( v );
}






