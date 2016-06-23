/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#ifndef BOGUS_EIGEN_PROBLEM_TRAITS_HPP
#define BOGUS_EIGEN_PROBLEM_TRAITS_HPP

#include "EigenMatrixTraits.hpp"

namespace bogus
{

template < typename LocalMatrixType >
struct ProblemTraits : public MatrixTraits< LocalMatrixType >
{
	typedef MatrixTraits< LocalMatrixType > Base ;
	typedef typename Base::Scalar Scalar ;
	enum{ dimension = Base::dimension } ;

	typedef Eigen::Matrix< Scalar, Eigen::Dynamic, 1 > DynVector ;
	typedef Eigen::Matrix< Scalar, Eigen::Dynamic, Eigen::Dynamic > DynMatrix ;
	template < typename OtherMatrix >
	struct MutableClone
	{
		typedef Eigen::Matrix< Scalar,
			OtherMatrix::RowsAtCompileTime,
			OtherMatrix::ColsAtCompileTime > Type ;
	} ;

	template< typename VectorType >
	static typename VectorType::template FixedSegmentReturnType< dimension >::Type segment( const unsigned i, VectorType& v )
	{ return v.template segment< dimension > ( i * dimension ) ; }
	template< typename VectorType >
	static typename VectorType::template ConstFixedSegmentReturnType< dimension >::Type segment( const unsigned i, const VectorType& v )
	{ return v.template segment< dimension > ( i * dimension ) ; }

} ;

template< unsigned Dimension, typename Scalar >
struct LocalProblemTraits : public ProblemTraits< Eigen::Matrix< Scalar, Dimension, Dimension > >
{

	typedef Eigen::Matrix< Scalar, Dimension, 1 > Vector ;
	typedef Eigen::Matrix< Scalar, Dimension, Dimension > Matrix ;

	typedef Eigen::Matrix< Scalar, Dimension-1, Dimension-1 > TgMatrix ;

	static Scalar np( const Vector & v )
	{ return v[0] ; }
	static Scalar& np( Vector & v )
	{ return v[0] ; }

	static typename Vector::template ConstFixedSegmentReturnType< Dimension - 1 >::Type
	tp( const Vector & v )  { return v.template segment< Dimension - 1 >( 1 ) ; }
	static typename Vector::template FixedSegmentReturnType< Dimension - 1 >::Type
	tp(       Vector & v )  { return v.template segment< Dimension - 1 >( 1 ) ; }

	static typename Matrix::ColXpr
	nc(      Matrix & m )  { return m.col( 0 ) ; }
	static Eigen::Block< Matrix, Dimension, Dimension -1 >
	tc(      Matrix & m )  { return m.template block< Dimension, Dimension - 1 >( 0, 1 ) ; }

	static Scalar&
	nnb(      Matrix & m )  { return m( 0, 0 ) ; }
	static Eigen::Block< Matrix, Dimension - 1, Dimension -1 >
	ttb(      Matrix & m )  { return m.template block< Dimension - 1, Dimension - 1 >( 1, 1 ) ; }
	static Eigen::Block< Matrix, Dimension - 1, 1 >
	tnb(      Matrix & m )  { return m.template block< Dimension - 1, 1 >( 1, 0 ) ; }
	static Eigen::Block< Matrix, 1, Dimension -1 >
	ntb(      Matrix & m )  { return m.template block< 1, Dimension - 1 >( 0, 1 ) ; }
} ;

template< typename Scalar >
struct LocalProblemTraits< Eigen::Dynamic, Scalar > : public ProblemTraits< Eigen::Matrix< Scalar, Eigen::Dynamic, Eigen::Dynamic > >
{
	typedef Eigen::Matrix< Scalar, Eigen::Dynamic, Eigen::Dynamic > Matrix ;
	typedef ProblemTraits< Matrix > Base ;

	typedef typename Base::DynVector Vector ;
	typedef typename Base::DynMatrix TgMatrix ;

	static Scalar np( const Vector & v )
	{ return v[0] ; }
	static Scalar& np( Vector & v )
	{ return v[0] ; }

	static typename Vector::ConstSegmentReturnType
	tp( const Vector & v )  { return v.segment( 1, v.rows() - 1 ) ; }
	static typename Vector::SegmentReturnType
	tp(       Vector & v )  { return v.segment( 1, v.rows() - 1 ) ; }

	static typename Matrix::ColXpr
	nc(      Matrix & m )  { return m.col( 0 ) ; }
	static Eigen::Block< Matrix >
	tc(      Matrix & m )  { return m.block( 0, 1, m.rows(), m.cols()-1 ) ; }

	static Scalar&
	nnb(      Matrix & m )  { return m( 0, 0 ) ; }
	static Eigen::Block< Matrix >
	ttb(      Matrix & m )  { return m.block( 1, 1, m.rows()-1, m.cols()-1 ) ; }
	static Eigen::Block< Matrix >
	tnb(      Matrix & m )  { return m.block( 1, 0, m.rows()-1, 1 ) ; }
	static Eigen::Block< Matrix >
	ntb(      Matrix & m )  { return m.block( 0, 1, 1, m.cols()-1 ) ; }
} ;


} //namespace bogus

#endif
