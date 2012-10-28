#ifndef SVMMATRIX
#define SVMMATRIX

#include "SparseGenMatrixT.h"

typedef SparseGenMatrixT<float> SvmMatrix;
typedef SmartPointer<SvmMatrix> SvmMatrixHandle;

#endif