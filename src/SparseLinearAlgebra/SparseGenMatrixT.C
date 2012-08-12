/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "SparseGenMatrixT.h"
#include "SparseStorageT.h"
#include <cassert>
#include "SimpleVector.h"

#include "DoubleMatrixTypes.h"

template <typename SCALAR>
int SparseGenMatrixT<SCALAR>::isKindOf( int type )
{
  return type == kSparseGenMatrix || type == kGenMatrix;
}


template<typename SCALAR>
SparseGenMatrixT<SCALAR>::SparseGenMatrixT( int rows, int cols, int nnz )
{
  mStorage = SparseStorageHandleT( new SparseStorageT<SCALAR>( rows, cols, nnz ) );
}


template<typename SCALAR>
SparseGenMatrixT<SCALAR>::SparseGenMatrixT( int rows, int cols, int nnz,
					  int krowM[], int jcolM[],
					  SCALAR M[] )
{
  mStorage = SparseStorageHandleT( new SparseStorageT<SCALAR>( rows, cols,
							 nnz, krowM,
							 jcolM, M ) );
}


template<typename SCALAR>
void SparseGenMatrixT<SCALAR>::atPutDense( int row, int col, double * A, int lda,
				      int rowExtent, int colExtent )
{
  mStorage->atPutDense( row, col, A, lda, rowExtent, colExtent );
}


template<typename SCALAR>
void SparseGenMatrixT<SCALAR>::fromGetDense( int row, int col, double * A, int lda,
					int rowExtent, int colExtent )
{
  mStorage->fromGetDense( row, col, A, lda, rowExtent, colExtent );
}
  

template<typename SCALAR>
void SparseGenMatrixT<SCALAR>::fromGetSpRow( int row, int col,
				    double A[], int lenA,
				    int jcolA[], int& nnz,
				    int colExtent, int& info )
{
  mStorage->fromGetSpRow( row, col, A, lenA, jcolA, nnz,
			  colExtent, info );
}


template<typename SCALAR>
void SparseGenMatrixT<SCALAR>::putSparseTriple( int irow[], int len,
					   int jcol[], double A[], 
					   int& info )
{
  mStorage->putSparseTriple( irow, len, jcol, A, info );
}


template<typename SCALAR>
void SparseGenMatrixT<SCALAR>::writeToStream(ostream& out) const
{
  mStorage->writeToStream( out );
}


template<typename SCALAR>
void SparseGenMatrixT<SCALAR>::randomize( double alpha, double beta, double * seed )
{
  mStorage->randomize( alpha, beta, seed );
}


template<typename SCALAR>
void SparseGenMatrixT<SCALAR>::getDiagonal( OoqpVector& vec )
{
  mStorage->getDiagonal( vec );
}


template<typename SCALAR>
void SparseGenMatrixT<SCALAR>::setToDiagonal( OoqpVector& vec )
{
  mStorage->setToDiagonal( vec );
}


template<typename SCALAR>
void SparseGenMatrixT<SCALAR>::atPutSpRow( int row, double A[],
				      int lenA, int jcolA[], int& info )
{
  mStorage->atPutSpRow( row, A, lenA, jcolA, info );
}


template<typename SCALAR>
int SparseGenMatrixT<SCALAR>::numberOfNonZeros()
{
  return mStorage->numberOfNonZeros();
}


template<typename SCALAR>
void SparseGenMatrixT<SCALAR>::symmetrize( int& info ) 
{
  mStorage->symmetrize( info );
}


template<typename SCALAR>
void SparseGenMatrixT<SCALAR>::getSize( int& m, int& n )
{
  m = mStorage->m;
  n = mStorage->n;
}


template<typename SCALAR>
void SparseGenMatrixT<SCALAR>::atPutSubmatrix( int destRow, int destCol,
					  DoubleMatrix& M,
					  int srcRow, int srcCol,
					  int rowExtent, int colExtent )
{
  int i, k;
  int info, nnz;

  int *    ja = new int[colExtent];
  double * a  = new double[colExtent];

  nnz = 0;
  for ( i = 0; i < rowExtent; i++ ) {
    M.fromGetSpRow( srcRow + i, srcCol, a, colExtent, ja,
		     nnz, colExtent, info );
    for( k = 0; k < nnz; k++ ) {
      ja[k] += (destCol - srcCol);
    }
    mStorage->atPutSpRow( destRow + i, a, nnz, ja, info );
  }

  delete [] ja;
  delete [] a;
}


template<typename SCALAR>
void SparseGenMatrixT<SCALAR>::mult ( double beta,  OoqpVector& y_in,
				 double alpha, OoqpVector& x_in )
{
  SimpleVector & x = dynamic_cast<SimpleVector &>(x_in);
  SimpleVector & y = dynamic_cast<SimpleVector &>(y_in);
  
  assert( x.n == mStorage->n && y.n == mStorage->m );

  double *xv = 0, *yv = 0;

  if( x.n > 0 ) xv = &x[0];
  if( y.n > 0 ) yv = &y[0];

  mStorage->mult( beta, yv, 1, alpha, xv, 1 );
}


template<typename SCALAR>
void SparseGenMatrixT<SCALAR>::transMult ( double beta,   OoqpVector& y_in,
				      double alpha,  OoqpVector& x_in )
{
  SimpleVector & x = dynamic_cast<SimpleVector &>(x_in);
  SimpleVector & y = dynamic_cast<SimpleVector &>(y_in);
  
  assert( x.n == mStorage->m && y.n == mStorage->n );

  double *xv = 0, *yv = 0;

  if( x.n > 0 ) xv = &x[0];
  if( y.n > 0 ) yv = &y[0];

  mStorage->transMult( beta, yv, 1, alpha, xv, 1 );
}




template<typename SCALAR>
double SparseGenMatrixT<SCALAR>::abmaxnorm()
{
  return mStorage->abmaxnorm();
}


template<typename SCALAR>
void SparseGenMatrixT<SCALAR>::atPutDiagonal( int idiag, OoqpVector& vvec )
{
  SimpleVector & v = dynamic_cast<SimpleVector &>(vvec);

  mStorage->atPutDiagonal( idiag, v.elements(), 1, v.n );
}


template<typename SCALAR>
void SparseGenMatrixT<SCALAR>::fromGetDiagonal( int idiag, OoqpVector& vvec )
{
  mStorage->fromGetDiagonal( idiag, vvec );
}

template<typename SCALAR>
void SparseGenMatrixT<SCALAR>::ColumnScale( OoqpVector& vec )
{
  mStorage->ColumnScale( vec );
}

template<typename SCALAR>
void SparseGenMatrixT<SCALAR>::SymmetricScale( OoqpVector& vec )
{
  mStorage->SymmetricScale( vec );
}

template<typename SCALAR>
void SparseGenMatrixT<SCALAR>::RowScale( OoqpVector& vec )
{
  mStorage->ColumnScale( vec );
}

template<typename SCALAR>
void SparseGenMatrixT<SCALAR>::scalarMult( double num )
{
  mStorage->scalarMult( num );
}
