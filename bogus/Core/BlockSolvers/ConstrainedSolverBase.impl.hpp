/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef BOGUS_CONSTRAINED_SOLVER_BASE_IMPL_HPP
#define BOGUS_CONSTRAINED_SOLVER_BASE_IMPL_HPP

#include "ConstrainedSolverBase.hpp"

#include "../Block/Access.hpp"
#include "../Block/BlockMatrixBase.hpp"

namespace bogus {


template < template <typename> class SolverType, typename BlockMatrixType >
template < typename NSLaw, typename RhsT, typename ResT >
typename ConstrainedSolverBase< SolverType,BlockMatrixType >::Scalar
ConstrainedSolverBase< SolverType,BlockMatrixType >::eval( const NSLaw &law,
							const ResT &y, const RhsT &x ) const
{
	const Segmenter< NSLaw::dimension, const RhsT, typename BlockMatrixType::Index >
			xSegmenter( x, m_matrix->rowOffsets() ) ;
	const Segmenter< NSLaw::dimension, const ResT, typename BlockMatrixType::Index >
			ySegmenter( y, m_matrix->rowOffsets() ) ;

	typedef typename BlockMatrixTraits< BlockMatrixType >::Index Index ;

	const Index n = m_matrix->rowsOfBlocks() ;

	Scalar err = 0., lres ;
	typename NSLaw::Traits::Vector lx, ly ;

	if( m_useInfinityNorm )
	{

	#ifndef BOGUS_DONT_PARALLELIZE
	#pragma omp parallel private( lx, ly, lres )
	#endif
		{
			lres = 0. ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp for
#endif
			for( Index i = 0 ; i < n ; ++ i )
			{
				lx = xSegmenter[ i ] * m_scaling[i] ;
				ly = ySegmenter[ i ] ;
				lres = std::max( law.eval( i, lx, ly ), lres ) ;
			}

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp critical
#endif
			err = std::max( err, lres ) ;
		}

		return err ;

	} else {

	#ifndef BOGUS_DONT_PARALLELIZE
	#pragma omp parallel for private( lx, ly, lres ) reduction ( + : err )
	#endif
		for( Index i = 0 ; i < n ; ++ i )
		{
			lx = xSegmenter[ i ] * m_scaling[i] ;
			ly = ySegmenter[ i ] ;
			// std::cout << "ly passed in as y: " << ly << std::endl;
//			std::cout << "ySegmenter workin with : " << ySegmenter << std::endl;
			lres = law.eval( i, lx, ly ) ;
			// std::cout << "ConstrainedSolverBaseImpl:: output lres " <<  lres << ",\t lx: " << lx << ",\t ly: " << ly << std::endl;
			err += lres ;
		}

		return err / ( 1 + n );

	}
}

template < template <typename> class SolverType, typename BlockMatrixType >
void ConstrainedSolverBase< SolverType,BlockMatrixType >::updateScalings()
{
	if( !m_matrix )
	{
		return ;
	}

	const Index n = m_matrix->rowsOfBlocks() ;
	m_scaling.resize( n ) ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
	for( Index i = 0 ; i <  n ; ++i )
	{
		m_scaling[i] = std::max( 1., m_matrix->diagonal(i).trace() ) ;
	}

}


} // namespace bogus

#endif
