/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#ifndef BOGUS_EIGEN_SPARSE_LINEAR_SOLVERS_HPP
#define BOGUS_EIGEN_SPARSE_LINEAR_SOLVERS_HPP

/*  Depending on your version of Eigen, this file may make use of LGPL licensed code.
	See Eigen/src/Sparse/SimplicialCholesky.h for more information
*/

#include "SparseHeader.hpp"

#ifdef BOGUS_WITH_EIGEN_STABLE_SPARSE_API

#include "../Utils/LinearSolverBase.hpp"
#include "../Utils/NaiveSharedPtr.hpp"

#include <Eigen/OrderingMethods>


#ifndef EIGEN_MPL2_ONLY

#include <Eigen/SparseCholesky>
#define BOGUS_WITH_EIGEN_SPARSE_LDLT

namespace bogus {

template < typename Derived >
struct LinearSolverTraits< LDLT< Eigen::SparseMatrixBase< Derived > > >
{
  typedef typename Derived::PlainObject MatrixType ;
  typedef Eigen::SimplicialLDLT< MatrixType > FactType ;

  template < typename RhsT > struct Result {
	  typedef Eigen::internal::solve_retval< Eigen::SimplicialCholeskyBase< FactType >, RhsT > Type ;
  } ;
  template < typename RhsT >
  struct Result< Eigen::MatrixBase< RhsT > > {
	  typedef typename Result< RhsT >::Type Type ;
  } ;
} ;


template < typename Derived >
struct LDLT< Eigen::SparseMatrixBase< Derived > >
		: public LinearSolverBase< LDLT< Eigen::SparseMatrixBase< Derived > > >
{
	typedef Eigen::SparseMatrixBase< Derived > MatrixType ;
	typedef LinearSolverTraits< LDLT< MatrixType > > Traits ;

	LDLT() {}
	template< typename OtherDerived >
	explicit LDLT ( const Eigen::SparseMatrixBase< OtherDerived >& mat )
		: m_fact( new typename Traits::FactType( mat ) )
	{}

	template< typename OtherDerived >
	LDLT& compute ( const Eigen::SparseMatrixBase< OtherDerived >& mat )
	{
		m_fact.reset( new typename Traits::FactType( mat ) ) ;
		return *this ;
	}

	template < typename RhsT, typename ResT >
	void solve( const Eigen::MatrixBase< RhsT >& rhs, ResT& res ) const
	{
		assert( m_fact ) ;
		res = m_fact->solve( rhs ) ;
	}

	template < typename RhsT >
	typename Traits::template Result< Eigen::MatrixBase< RhsT > >::Type
	solve( const Eigen::MatrixBase< RhsT >& rhs ) const
	{
		assert( m_fact ) ;
		return m_fact->solve( rhs ) ;
	}

  private:
	BOGUS_SHARED_PTR( typename Traits::FactType, m_fact ) ;
} ;

template < typename Scalar, int _Options = 0, typename _Index = int >
struct SparseLDLT : public LDLT< Eigen::SparseMatrixBase< Eigen::SparseMatrix< Scalar, _Options, _Index > > >
{
	SparseLDLT() {}
	template< typename OtherDerived >
	explicit SparseLDLT ( const Eigen::SparseMatrixBase< OtherDerived >& mat )
		: LDLT< Eigen::SparseMatrixBase< Eigen::SparseMatrix< Scalar, _Options, _Index > > >( mat )
	{}
} ;

} //snamespace bogus

#endif // LDLT


#if EIGEN_VERSION_AT_LEAST(3,1,92)

#include <Eigen/SparseLU>

#define BOGUS_WITH_EIGEN_SPARSE_LU

namespace bogus {

template < typename Derived >
struct LinearSolverTraits< LU< Eigen::SparseMatrixBase< Derived > > >
{
  typedef typename Derived::PlainObject MatrixType ;
  typedef Eigen::SparseLU< MatrixType, Eigen::COLAMDOrdering<int> > FactType ;

  template < typename RhsT > struct Result {
	  typedef Eigen::internal::solve_retval< FactType, RhsT > Type ;
  } ;
  template < typename RhsT >
  struct Result< Eigen::MatrixBase< RhsT > > {
	  typedef typename Result< RhsT >::Type Type ;
  } ;
} ;


template < typename Derived >
struct LU< Eigen::SparseMatrixBase< Derived > >
		: public LinearSolverBase< LU< Eigen::SparseMatrixBase< Derived > > >
{
	typedef Eigen::SparseMatrixBase< Derived > MatrixType ;
	typedef LinearSolverTraits< LU< MatrixType > > Traits ;

	LU() {}
	template< typename OtherDerived >
	explicit LU ( const Eigen::SparseMatrixBase< OtherDerived >& mat )
		: m_fact( new typename Traits::FactType( mat ) )
	{}

	template< typename OtherDerived >
	LU& compute ( const Eigen::SparseMatrixBase< OtherDerived >& mat )
	{
		m_fact.reset( new typename Traits::FactType( mat ) ) ;
		return *this ;
	}

	template < typename RhsT, typename ResT >
	void solve( const Eigen::MatrixBase< RhsT >& rhs, ResT& res ) const
	{
		assert( m_fact ) ;
		res = m_fact->solve( rhs ) ;
	}

	template < typename RhsT >
	typename Traits::template Result< Eigen::MatrixBase< RhsT > >::Type
	solve( const Eigen::MatrixBase< RhsT >& rhs ) const
	{
		assert( m_fact ) ;
		return m_fact->solve( rhs ) ;
	}

  private:
	BOGUS_SHARED_PTR( typename Traits::FactType, m_fact ) ;
} ;

template < typename Scalar, int _Options = 0, typename _Index = int >
struct SparseLU : public LU< Eigen::SparseMatrixBase< Eigen::SparseMatrix< Scalar, _Options, _Index > > >
{
	SparseLU() {}
	template< typename OtherDerived >
	explicit SparseLU ( const Eigen::SparseMatrixBase< OtherDerived >& mat )
		: LU< Eigen::SparseMatrixBase< Eigen::SparseMatrix< Scalar, _Options, _Index > > >( mat )
	{}
} ;

} //namespace bogus

#endif // LU

#endif // EIGEN_STABLE_API

#endif //HPP
