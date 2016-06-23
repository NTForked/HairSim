/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#ifndef BOGUS_BLOCK_GAUSS_SEIDEL_BASE_IMPL_HPP
#define BOGUS_BLOCK_GAUSS_SEIDEL_BASE_IMPL_HPP

#include "GaussSeidelBase.hpp"
#include "ConstrainedSolverBase.impl.hpp"

#ifndef BOGUS_DONT_PARALLELIZE
#include <omp.h>
#endif

namespace bogus
{

template < template <typename> class GaussSeidelImpl, typename BlockMatrixType >
void GaussSeidelBase< GaussSeidelImpl, BlockMatrixType >::updateLocalMatrices( )
{

	Base::updateScalings() ;

	if( !m_matrix )
		return ;

	const Index n = m_matrix->rowsOfBlocks() ;
	m_localMatrices.resize( n ) ;
	m_regularization.resize( n ) ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
	for( Index i = 0 ; i <  n ; ++i )
	{
		m_localMatrices[i] = m_matrix->diagonal( i ) ;

		if( m_autoRegularization > 0. )
		{
			m_regularization(i) = std::max( 0., m_autoRegularization - m_localMatrices[i].eigenvalues().real().minCoeff() ) ;
			m_localMatrices[i].diagonal().array() += m_regularization(i) ;
		} else m_regularization(i) = 0. ;

		// std::cout << "m_localMatrices[" << i << "] is a [ " << m_localMatrices[i].rows() << " x " << m_localMatrices[i].cols() << " ] matrix: " << std::endl
		//  << m_localMatrices[i](0,0) << ", " << m_localMatrices[i](0,1) << ", " << m_localMatrices[i](0,2) << ";" << std::endl 
		//  << m_localMatrices[i](1,0) << ", " << m_localMatrices[i](1,1) << ", " << m_localMatrices[i](1,2) << ";" << std::endl 
		//  << m_localMatrices[i](2,0) << ", " << m_localMatrices[i](2,1) << ", " << m_localMatrices[i](2,2) << ";" << std::endl; 

	}
}


} //namespace bogus


#endif
