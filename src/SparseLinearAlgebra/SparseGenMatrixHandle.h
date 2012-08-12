/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef SPARSEGENMATRIXHANDLE_H
#define  SPARSEGENMATRIXHANDLE_H

#include "IotrRefCount.h"
#include "SmartPointer.h"

template <typename S> class SparseGenMatrixT;
typedef SparseGenMatrixT<double> SparseGenMatrix;
typedef SmartPointer<SparseGenMatrixT<double> > SparseGenMatrixHandle;


#endif
