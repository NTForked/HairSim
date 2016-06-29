#ifndef BANDMATRIXFWD_H
#define BANDMATRIXFWD_H

#include "../Utils/Definitions.h"

template<typename ScalarT, IndexType kl, IndexType ku>
class BandMatrix;

typedef BandMatrix<Scalar, 10, 10> JacobianMatrixType;
typedef BandMatrix<Scalar, 15, 15> RestJacobianMatrixType;
typedef BandMatrix<Scalar, 1, 1> TriDiagonalMatrixType;

template<typename ScalarT, int kl > 
class SymmetricBandMatrixSolver;

typedef SymmetricBandMatrixSolver<Scalar, 10> JacobianSolver;

#endif // BANDMATRIXFWD_HH
