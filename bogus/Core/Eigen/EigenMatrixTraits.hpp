/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#ifndef BOGUS_EIGEN_MATRIX_TRAITS_HPP
#define BOGUS_EIGEN_MATRIX_TRAITS_HPP

#include "../Utils/LinearSolverBase.hpp"

#include <Eigen/Core>
#include "SparseHeader.hpp"

namespace bogus
{

template < typename _MatrixType >
struct MatrixTraits
{
	typedef _MatrixType MatrixType ;
	enum { dimension = MatrixType::RowsAtCompileTime  } ;
	typedef typename MatrixType::Scalar Scalar ;

	typedef LU  < Eigen::MatrixBase< MatrixType > > LUType ;
	typedef LDLT< Eigen::MatrixBase< MatrixType > > LDLTType ;

} ;

template < typename _Scalar, int _Options, typename _Index >
struct MatrixTraits< Eigen::SparseMatrix< _Scalar, _Options, _Index > >
{
	typedef _Scalar Scalar ;
	enum { dimension = -1  } ;
	typedef Eigen::SparseMatrix< Scalar, _Options, _Index > MatrixType ;

	typedef LU< Eigen::SparseMatrixBase< Eigen::SparseMatrix< Scalar, _Options, _Index > > > LUType ;
	typedef LDLT< Eigen::SparseMatrixBase< Eigen::SparseMatrix< Scalar, _Options, _Index > > > LDLTType ;

} ;

}

#endif
