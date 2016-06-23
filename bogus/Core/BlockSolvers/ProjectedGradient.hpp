/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#ifndef BOGUS_BLOCK_PROJECTED_GRADIENT_HPP
#define BOGUS_BLOCK_PROJECTED_GRADIENT_HPP

#include "ConstrainedSolverBase.hpp"

#include <vector>

namespace bogus
{

//! Projected Gradient iterative solver.
template < typename BlockMatrixType >
class ProjectedGradient : public ConstrainedSolverBase< ProjectedGradient, BlockMatrixType >
{
public:
	typedef ConstrainedSolverBase< bogus::ProjectedGradient, BlockMatrixType > Base ;

	typedef typename Base::GlobalProblemTraits GlobalProblemTraits ;
	typedef typename GlobalProblemTraits::Scalar Scalar ;

	//! Default constructor -- you will have to call setMatrix() before using the solve() function
	ProjectedGradient( ) : Base() { init() ; }
	//! Constructor with the system matrix
	explicit ProjectedGradient( const BlockMatrixBase< BlockMatrixType > & matrix ) : Base()
	{ init() ; Base::setMatrix( matrix ) ; }

	//! Finds an approximate solution for a constrained linear problem
	template < typename NSLaw, typename RhsT, typename ResT >
	Scalar solve( const NSLaw &law, const RhsT &b, ResT &x ) const ;


	void setMatrix( const BlockMatrixBase< BlockMatrixType > & matrix )
	{
		m_matrix = &matrix ;
		Base::updateScalings() ;
	}

protected:

	typedef typename Base::Index Index ;

	template < typename NSLaw, typename VectorT >
	void projectOnConstraints( const NSLaw &projector, VectorT &x ) const ;

	//! Sets up the default values for all parameters
	void init()
	{
		m_tol = 1.e-6 ;
		m_maxIters = 300 ;
		m_lsIters = 8 ;
		m_lsOptimisticFactor = 2 ;
		m_lsPessimisticFactor = .5 ;
		m_lsArmijoCriterion = 1.e-4;
	}

	using Base::m_matrix ;
	using Base::m_maxIters ;
	using Base::m_tol ;

	Scalar m_lsIters ;
	Scalar m_lsOptimisticFactor ;
	Scalar m_lsPessimisticFactor ;
	Scalar m_lsArmijoCriterion ;

} ;

} //namespace bogus


#endif

