/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#ifndef BOGUS_SPARSEBLOCKMATRIX_HPP
#define BOGUS_SPARSEBLOCKMATRIX_HPP

#include "SparseBlockMatrixBase.hpp"

namespace bogus
{

//! Specialization of BlockMatrixTraits for SparseBlockMatrix
template < typename BlockT, int Flags >
struct BlockMatrixTraits< SparseBlockMatrix< BlockT, Flags > >
		: public BlockMatrixTraits< BlockObjectBase< SparseBlockMatrix< BlockT, Flags > > >
{
	typedef BlockMatrixTraits< BlockObjectBase< SparseBlockMatrix< BlockT, Flags > > > BaseTraits ;
	typedef typename BaseTraits::Index      Index;
	typedef typename BaseTraits::BlockPtr   BlockPtr;

	typedef BlockT BlockType ;
	typedef typename BlockTraits< BlockT >::Scalar Scalar ;
	typedef typename ResizableSequenceContainer< BlockType >::Type BlocksArrayType ;

	enum {
		is_transposed = 0,
		is_temporary  = 0,

		is_compressed = ~Flags & flags::UNCOMPRESSED,
		is_symmetric  =  Flags & flags::SYMMETRIC,
		is_col_major  =  Flags & flags::COL_MAJOR,
		flags         =  Flags
	} ;

	typedef SparseBlockIndex< is_compressed, Index, BlockPtr > MajorIndexType ;

} ;

//! Sparse Block Matrix
/*!
  \tparam BlockT the type of the blocks of the matrix. Can be scalar, Eigen dense of sparse matrices,
  or basically anything provided a few functions are specialized
  \tparam Flags a combination of the values defined in \ref bogus::flags
  */
template < typename BlockT, int Flags >
class SparseBlockMatrix : public  SparseBlockMatrixBase< SparseBlockMatrix< BlockT, Flags > >
{
public:
	typedef SparseBlockMatrixBase< SparseBlockMatrix > Base ;

	SparseBlockMatrix() : Base()
	{
	}

	template < typename RhsT >
	SparseBlockMatrix( const BlockObjectBase< RhsT >& rhs ) : Base()
	{
		Base::operator= ( rhs.derived() ) ;
	}

	template < typename RhsT >
	SparseBlockMatrix& operator=( const BlockObjectBase< RhsT >& rhs )
	{
		return ( Base::operator= ( rhs.derived() ) ).derived() ;
	}

} ;

// Specialization for block matrix of MappedSparseBlockMatrix
template < typename BlockT, int Flags >
struct BlockTraits< SparseBlockMatrix< BlockT, Flags > >
{
	typedef SparseBlockMatrix< BlockT, Flags > BlockType ;

	typedef typename BlockType::Scalar Scalar ;
	typedef SparseBlockMatrix
		< typename BlockType::TransposeBlockType, Flags ^ flags::COL_MAJOR >
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
