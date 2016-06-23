/*
 * This file is part of So-bogus, a C++ sparse block matrix library and
 * Second Order Cone solver.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * So-bogus is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.

 * So-bogus is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with So-bogus.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "FrictionProblem.hpp"

#include "../Core/Block.impl.hpp"
#include "../Extra/SecondOrder.impl.hpp"

#include "../Extra/SecondOrder.impl.hpp"

namespace bogus {

namespace friction_problem {

template< unsigned Dimension, typename EigenDerived, typename Index >
void applyPermutation(
		const std::vector< std::size_t >& permutation,
		Eigen::MatrixBase< EigenDerived > &vec,
		const Index* offsets
		)
{
	Segmenter< Dimension, EigenDerived, Index > segmenter( vec.derived(), offsets ) ;
	bogus::applyPermutation( permutation.size(), &permutation[0], segmenter ) ;
}


template< unsigned Dimension, template <typename> class Method >
static double solve( const DualFrictionProblem< Dimension >& dual,
		const ConstrainedSolverBase< Method, typename DualFrictionProblem< Dimension >::WType > &gs,
		double *r, const bool staticProblem )
{
	typename Eigen::VectorXd::MapType r_map ( r, dual.W.rows() ) ;

	if( dual.permuted() )
		applyPermutation< Dimension >( dual.permutation(), r_map, dual.W.majorIndex().innerOffsetsData() ) ;

	double res = staticProblem
			? gs.solve( typename DualFrictionProblem< Dimension >::SOCLawType
						( dual.W.rowsOfBlocks(), dual.mu.data() ), dual.b, r_map )
			: gs.solve( typename DualFrictionProblem< Dimension >::CoulombLawType
						( dual.W.rowsOfBlocks(), dual.mu.data() ), dual.b, r_map ) ;

	if( dual.permuted() )
		applyPermutation< Dimension >( dual.invPermutation(), r_map, dual.W.majorIndex().innerOffsetsData() ) ;

	return res ;
}

template< unsigned Dimension, template <typename> class Method >
static double eval( const DualFrictionProblem< Dimension >& dual,
		const ConstrainedSolverBase< Method, typename DualFrictionProblem< Dimension >::WType > &gs,
		const double *r_data, const bool staticProblem  )
{
	Eigen::VectorXd r = Eigen::VectorXd::Map( r_data, dual.W.rows() ) ;
	Eigen::VectorXd u = dual.W*r + dual.b ;

	if( dual.permuted())
		applyPermutation< Dimension >( dual.permutation(), r, dual.W.majorIndex().innerOffsetsData() ) ;


	double res = staticProblem
			? gs.eval( typename DualFrictionProblem< Dimension >::SOCLawType
					   ( dual.W.rowsOfBlocks(), dual.mu.data() ), u, r )
			: gs.eval( typename DualFrictionProblem< Dimension >::CoulombLawType
					   ( dual.W.rowsOfBlocks(), dual.mu.data() ), u, r ) ;

	return res ;
}

template< unsigned Dimension, template <typename> class Method >
static double solveCadoux( const DualFrictionProblem< Dimension >& dual,
		const ConstrainedSolverBase< Method, typename DualFrictionProblem< Dimension >::WType > &gs,
		double *r, const unsigned cadouxIterations, const Signal<unsigned, double> *callback )
{
	const std::ptrdiff_t n = dual.W.rowsOfBlocks() ;

	typename DualFrictionProblem< Dimension >::CoulombLawType coulombLaw( n, dual.mu.data() ) ;
	typename DualFrictionProblem< Dimension >::SOCLawType         socLaw( n, dual.mu.data() ) ;

	gs.setMatrix( dual.W );
	Eigen::Map< Eigen::VectorXd > r_map ( r, dual.W.rows() ) ;

	if( dual.permuted() )
		applyPermutation< Dimension >( dual.permutation(), r_map, dual.W.majorIndex().innerOffsetsData() ) ;

	Eigen::VectorXd s( dual.W.rows() ) ;

	double res = -1 ;
	const double tol = gs.tol() ;
	gs.setTol( 1.e-1 * tol ) ;	//dual.We might experience slow convergence is GS not precise enough

	for( unsigned cdxIter = 0 ; cdxIter < cadouxIterations ; ++cdxIter )
	{
		s = dual.W * r_map + dual.b ;

		res = gs.eval( coulombLaw, s, r_map ) ;

		if( callback ) callback->trigger( cdxIter, res ) ;
		if( cdxIter > 0 && res < tol ) break ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
		for( std::ptrdiff_t i = 0 ; i < n ; ++i )
		{
			s[ Dimension*i ] = s.segment< Dimension-1 >( Dimension*i+1 ).norm() * dual.mu[i] ;
			s.segment< Dimension-1  >( Dimension*i+1 ).setZero() ;
		}

		s += dual.b ;

		gs.solve( socLaw, s, r_map ) ;

	}

	gs.setTol( tol ) ;

	if( dual.permuted() )
		applyPermutation< Dimension >( dual.invPermutation(), r_map, dual.W.majorIndex().innerOffsetsData() ) ;

	return res ;
}

} //namespace friction_problem

} //namespace bogus
