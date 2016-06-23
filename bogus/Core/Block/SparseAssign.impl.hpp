/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#ifndef BOGUS_SPARSEASSIGN_IMPL_HPP
#define BOGUS_SPARSEASSIGN_IMPL_HPP

#include "Access.hpp"

#include "SparseBlockMatrix.hpp"
#include "SparseBlockIndexComputer.hpp"

namespace bogus
{

template< bool Transpose >
struct BlockCopier
{
	template < typename BlockT1, typename BlockT2, typename ScalarT >
	static void copy( const BlockT2* source, BlockT1* dest, int n, ScalarT scale )
	{
		typedef TransposeIf< Transpose > Getter ;
		if( scale == 1 )
		{
#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
			for( int i = 0 ; i < n ; ++i )
				dest[i] = Getter::get( source[i] ) ;
		} else {
#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
			for( int i = 0 ; i < n ; ++i )
				dest[i] = scale * Getter::get( source[i] ) ;
		}
	}
} ;

template < typename Derived >
Derived& SparseBlockMatrixBase<Derived>::operator= ( const SparseBlockMatrixBase &source )
{
	m_rows = source.rows() ;
	m_cols = source.cols() ;

	m_blocks = source.blocks() ;
	m_transposeBlocks = source.transposeBlocks() ;

	m_majorIndex = source.majorIndex() ;
	m_minorIndex = source.minorIndex() ;
	m_transposeIndex = source.transposeIndex() ;

	return derived() ;
}

template < typename Derived >
template < bool Transpose, typename OtherDerived >
Derived& SparseBlockMatrixBase<Derived>::assign( const SparseBlockMatrixBase< OtherDerived > &source, Scalar scale )
{
	BOGUS_STATIC_ASSERT( !Transpose || IsTransposable< typename OtherDerived::BlockType >::Value,
						 TRANSPOSE_IS_NOT_DEFINED_FOR_THIS_BLOCK_TYPE
	) ;

	if( static_cast< const void* >( this ) == static_cast< const void* >( &source ) ) return derived() ;

	typedef typename SparseBlockMatrixBase< OtherDerived >::Traits OtherTraits ;
	const bool sameMajorness = Transpose ^ ( ( (bool) Traits::is_col_major )  == ( (bool) OtherTraits::is_col_major ) ) ;
	const bool sameSymmetry =  ( (bool) Traits::is_symmetric )  == ( (bool) OtherTraits::is_symmetric ) ;
	bool useTransposeIndex = false ;

	if( ( Traits::is_symmetric && OtherTraits::is_symmetric ) || sameMajorness )
	{
		m_majorIndex = source.majorIndex() ;
		m_minorIndex = source.minorIndex() ;
		m_transposeIndex = source.transposeIndex() ;
	}
	else if ( Transpose && source.transposeIndex().valid )
	{
		m_majorIndex = source.transposeIndex() ;
		m_transposeIndex = source.majorIndex() ;
		m_minorIndex.innerOffsets = m_transposeIndex.innerOffsets ;
		m_minorIndex.clear() ;
		m_minorIndex.valid = false ;
		useTransposeIndex = true ;
	} else {
		if( Transpose )
		{
			rowMajorIndex() = source.colMajorIndex() ;
			colMajorIndex() = source.rowMajorIndex() ;
			m_transposeIndex.valid = false ;
		} else {
			rowMajorIndex() = source.rowMajorIndex() ;
			colMajorIndex() = source.colMajorIndex() ;
			m_transposeIndex = source.transposeIndex() ;
		}
	}

	m_majorIndex.valid &= sameSymmetry ;

	m_cols = source.cols() ;
	m_rows = source.rows() ;
	if( Transpose ) std::swap( m_cols, m_rows ) ;

	if( m_majorIndex.valid )
	{

		if( useTransposeIndex )
		{
			assert( m_transposeIndex.valid ) ;
			m_blocks.resize( source.transposeBlocks().size() ) ;
			m_transposeBlocks.resize( source.blocks().size() ) ;
			BlockCopier< false >::copy( source.transposeData(), this->data(), source.transposeBlocks().size(), scale ) ;
			BlockCopier< false >::copy( source.data(), this->transposeData(), source.blocks().size(), scale ) ;
		} else {
			const bool needTranspose = Transpose != ( Traits::is_symmetric && OtherTraits::is_symmetric && !sameMajorness ) ;

			m_blocks.resize( source.blocks().size() ) ;
			BlockCopier< needTranspose >::copy( source.data(), this->data(), source.blocks().size(), scale ) ;

			if( m_transposeIndex.valid )
			{
				m_transposeBlocks.resize( source.transposeBlocks().size() ) ;
				BlockCopier< needTranspose >::copy( source.transposeData(), this->transposeData(), source.transposeBlocks().size(), scale ) ;
			}
		}

		Finalizer::finalize( *this ) ;
	} else {

		clear() ;
		reserve( source.blocks().size() ) ;

		typedef SparseBlockIndexComputer< OtherDerived, Traits::is_col_major, Transpose > IndexComputerType ;
		IndexComputerType indexComputer( source ) ;
		typedef typename IndexComputerType::ReturnType SourceIndexType ;
		const SourceIndexType &sourceIndex = indexComputer.get() ;

		typedef BlockTransposeOption<
				OtherTraits::is_symmetric && !( BlockTraits< typename OtherTraits::BlockType >::is_self_transpose ),
				Transpose > TransposeIf ;

		assert( sourceIndex.valid ) ;

		for( Index i = 0 ; i < sourceIndex.outerSize() ; ++i )
		{
			for( typename SourceIndexType::InnerIterator src_it( sourceIndex, i ) ;
				 src_it && !( Traits::is_symmetric && i < src_it.inner() ) ; ++ src_it )
			{
				TransposeIf::assign( source.block( src_it.ptr() ) ,
									 insertByOuterInner< true >( i, src_it.inner() ),
									 scale, src_it.after( i ) ) ;
			}
		}
		finalize() ;
	}


	return derived();
}

} //namespace bogus

#endif // SPARSETRANSPOSE_IMPL_HPP
