/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef BOGUS_PROJECTED_GRADIENT_IMPL_HPP
#define BOGUS_PROJECTED_GRADIENT_IMPL_HPP


#include "ProjectedGradient.hpp"
#include "ConstrainedSolverBase.impl.hpp"

namespace bogus
{

template < typename BlockMatrixType >
template < typename NSLaw, typename RhsT, typename ResT >
typename ProjectedGradient< BlockMatrixType >::Scalar
ProjectedGradient< BlockMatrixType >::solve(
		const NSLaw &law, const RhsT &b, ResT &x ) const
{

	typename GlobalProblemTraits::DynVector
			Mx ( b.rows() ),
			y  ( b.rows() ), // = Mx +b
			xs ( x.rows() )  // tentative new value for x
			;

	projectOnConstraints( law, x ) ;

	// Unconstrained objective function
	Mx = (*m_matrix)*x ;
	Scalar J = x.dot( .5 * Mx + b ) ;

	Scalar res = -1, alpha = 1 ;

	for( unsigned pgIter = 0 ; pgIter < m_maxIters ; ++pgIter )
	{
		// y = grad J = Mx+b
		y = Mx + b ;
		res = Base::eval( law, y, x ) ;

		this->m_callback.trigger( pgIter, res );
		if( res < m_tol ) break ;

		Scalar Js = J ;

		alpha *= m_lsOptimisticFactor ;

		// Armijo line-search
		for( unsigned lsIter = 0 ;
			 lsIter < m_lsIters ;
			 ++ lsIter, alpha *= m_lsPessimisticFactor )
		{
			xs = x - alpha * y ;
			projectOnConstraints( law, xs ) ;

			Mx = (*m_matrix) * xs ;
			Js = xs.dot( .5 * Mx + b ) ;

			if( Js < J + m_lsArmijoCriterion * y.dot( xs - x ) )
				break ;

		}

		x = xs ;
		J = Js ;
	}


	return res ;
}

template < typename BlockMatrixType >
template < typename NSLaw, typename VectorT >
void ProjectedGradient< BlockMatrixType >::projectOnConstraints(
		const NSLaw &law, VectorT &x ) const
{
	Segmenter< NSLaw::dimension, VectorT, typename BlockMatrixType::Index >
			xSegmenter( x, m_matrix->rowOffsets() ) ;

	const Index n = m_matrix->rowsOfBlocks() ;
	typename NSLaw::Traits::Vector lx ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for private( lx )
#endif
	for( Index i = 0 ; i < n ; ++ i )
	{
		lx = xSegmenter[ i ] ;
		law.projectOnConstraint( i, lx ) ;
		xSegmenter[ i ] = lx ;
	}

}

} //namespace bogus

#endif
