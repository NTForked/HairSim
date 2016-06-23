/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#ifndef BOGUS_BLOCKOBJECTBASE_HPP
#define BOGUS_BLOCKOBJECTBASE_HPP

#include "../Block.fwd.hpp"

namespace bogus
{

//! Base class for anything block
template < typename Derived >
struct BlockObjectBase
{
	//! Returns a const reference to the implementation
	const Derived& derived() const
	{
		return static_cast< const Derived& >( *this ) ;
	}
	//! Returns a reference to the implementation
	Derived& derived()
	{
		return static_cast< Derived& >( *this ) ;
	}

	typedef BlockMatrixTraits< Derived > Traits ;

	typedef typename Traits::Index Index ;
	typedef typename Traits::Scalar Scalar ;
	typedef typename Traits::ConstTransposeReturnType ConstTransposeReturnType ;
	typedef typename Traits::TransposeObjectType TransposeObjectType ;

	typedef typename Traits::PlainObjectType PlainObjectType ;
	typedef typename Traits::EvalType EvalType ;
	enum { is_transposed = Traits::is_transposed } ;

	//! Returns the total number of rows of the matrix ( expanding blocks )
	Index rows() const { return derived().rows() ; }
	//! Returns the total number of columns of the matrix ( expanding blocks )
	Index cols() const { return derived().cols() ; }

	//! Return a const transposed view of this object
	ConstTransposeReturnType transpose() const { return derived().transpose() ; }

	//! Eval this object in a temporary. For internal use, not part of the public API
	EvalType eval() const { return derived().eval() ; }
};

//! Default specialization of traits for BlockMatrices
/*! Re-specialiazed for derived classes, see e.g. BlockMatrixTraits< SparseBlockMatrix > */
template< typename Derived  >
struct BlockMatrixTraits< BlockObjectBase< Derived > > {

	typedef BOGUS_DEFAULT_INDEX_TYPE    Index ;
	typedef BOGUS_DEFAULT_BLOCK_PTR_TYPE BlockPtr ;

	typedef Derived PlainObjectType ;
	typedef const PlainObjectType* EvalType ;

	typedef Transpose< Derived > ConstTransposeReturnType ;
	typedef ConstTransposeReturnType TransposeObjectType ;
} ;

}

#endif // BLOCKMATRIX_HPP
