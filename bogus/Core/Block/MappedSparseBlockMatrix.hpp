/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#ifndef BOGUS_MAPPED_SPARSEBLOCKMATRIX_HPP
#define BOGUS_MAPPED_SPARSEBLOCKMATRIX_HPP

#include "CompressedSparseBlockIndex.hpp"
#include "SparseBlockMatrixBase.hpp"
#include "../Utils/CppTools.hpp"

namespace bogus
{

//! Specialization of BlockMatrixTraits for SparseBlockMatrix
template < typename BlockT, int Flags, typename Index_ >
struct BlockMatrixTraits< MappedSparseBlockMatrix< BlockT, Flags, Index_ > >
		: public BlockMatrixTraits< BlockObjectBase< MappedSparseBlockMatrix< BlockT, Flags, Index_ > > >
{
	typedef BlockMatrixTraits< BlockObjectBase< MappedSparseBlockMatrix< BlockT, Flags, Index_ > > > BaseTraits ;
	typedef Index_      Index;
	typedef typename BaseTraits::BlockPtr   BlockPtr;

	typedef BlockT BlockType ;
	typedef typename BlockTraits< BlockT >::Scalar Scalar ;
	typedef typename ConstMappedArray< BlockType >::Type BlocksArrayType ;

	enum {
		is_transposed  = 0,
		is_temporary   = 0,

		is_compressed  = 1,
		is_symmetric   = Flags & flags::SYMMETRIC,
		is_col_major   = Flags & flags::COL_MAJOR,
		flags          = Flags & ~flags::UNCOMPRESSED
	} ;

	typedef SparseBlockIndex< is_compressed, Index, BlockPtr, ConstMappedArray > MajorIndexType ;

} ;

//! Mapped Sparse Block Matrix
/*!

  Allows bogus to operate on an external, immutable matrix using a compressed index
  ( ala MKL BSR -- Block Sparse Row )
  Only const (or caching) operations are available.

  To obtain a valid MappedSparseBlockMatrix from a set of raw pointers, two operations are necessary
   1 Set the size of each row and column of blocks using the setRows()/setCols() or cloneDimensions() methods
   2 Map the block data and index using the mapTo() method

  \tparam BlockT the type of the blocks of the matrix. Can be scalar, Eigen dense of sparse matrices,
  or basically anything provided a few functions are specialized
  \tparam Flags a combination of the values defined in \ref bogus::flags
  \tparam The integer type used for indexing
  */
template < typename BlockT, int Flags, typename Index_ >
class MappedSparseBlockMatrix : public  SparseBlockMatrixBase< MappedSparseBlockMatrix< BlockT, Flags, Index_ > >
{
public:
	typedef SparseBlockMatrixBase< MappedSparseBlockMatrix > Base ;
	typedef typename Base::Index    Index ;
	typedef typename Base::BlockPtr BlockPtr ;

	MappedSparseBlockMatrix() : Base()
	{
	}

	template < typename Derived >
	explicit MappedSparseBlockMatrix( const SparseBlockMatrixBase< Derived > &source ) : Base()
	{
		mapTo( source ) ;
	}

	void mapTo(
			std::size_t numberOfNonZeros,    //!< Total number of blocks
			const BlockT* dataPtr,           //!< Pointer to an array containing the blocks data
			const Index*  outerIndexPtr,     //!< A.k.a rowsIndex, pntrb, pntre-1 from BSR format
			const Index*  innerIndexPtr      //!< A.k.a columns from BSR format
			)
	{
		Base::clear() ;

		Base::m_blocks.setData( dataPtr, numberOfNonZeros ) ;
		Base::m_majorIndex.inner.setData( innerIndexPtr, numberOfNonZeros ) ;
		Base::m_majorIndex.outer.setData( outerIndexPtr, Base::m_minorIndex.innerSize() + 1 ) ;

		Base::m_minorIndex.valid = Base::empty() ;
		Base::Finalizer::finalize( *this ) ;
	}

	template < typename Derived >
	void mapTo( const SparseBlockMatrixBase< Derived > &source )
	{
		BOGUS_STATIC_ASSERT( (int) Base::Traits::flags == (int) Derived::Traits::flags,
							 OPERANDS_HAVE_INCONSISTENT_FLAGS ) ;

		Base::cloneDimensions( source ) ;
		mapTo( source.nBlocks(),
			   source.data(),
			   source.majorIndex().outerIndexPtr(),
			   source.majorIndex().innerIndexPtr()
			   );

	}

} ;

// Specialization for block matrix of MappedSparseBlockMatrix
template < typename BlockT, int Flags, typename Index_ >
struct BlockTraits< MappedSparseBlockMatrix< BlockT, Flags, Index_ > >
{
	typedef MappedSparseBlockMatrix< BlockT, Flags, Index_ > BlockType ;
	typedef typename BlockType::Scalar Scalar ;
	typedef MappedSparseBlockMatrix
		< typename BlockType::TransposeBlockType, Flags ^ flags::COL_MAJOR, Index_ >
		TransposeStorageType ;

	enum {
		RowsAtCompileTime = internal::DYNAMIC,
		ColsAtCompileTime = internal::DYNAMIC,
		uses_plain_array_storage = 0,
		is_row_major = !BlockMatrixTraits< BlockType >::is_col_major,
		is_self_transpose = BlockMatrixTraits< BlockType >::is_symmetric
	}  ;
} ;

}

#endif // SPARSEBLOCKMATRIX_HH
