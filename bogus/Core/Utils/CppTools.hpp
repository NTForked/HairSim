/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef BOGUS_CPP_TOOLS_HPP
#define BOGUS_CPP_TOOLS_HPP

#include <iostream>

namespace bogus
{

// Swap template parameters if DoSwap is true

template < bool DoSwap, typename First_, typename Second_ >
struct TypeSwapIf
{
	typedef First_  First  ;
	typedef Second_ Second ;
} ;

template < typename First_, typename Second_ >
struct TypeSwapIf< true, First_, Second_ >
{
	typedef First_  Second  ;
	typedef Second_ First ;
} ;

template < bool DoSwap, int First_, int Second_ >
struct SwapIf
{
	enum { First = First_, Second = Second_  } ;
} ;

template < int First_, int Second_ >
struct SwapIf< true, First_, Second_ >
{
	enum { First = Second_, Second = First_  } ;
} ;


// Enable if (for SFINAE use )

template < bool Condition, typename ReturnType_ = void >
struct EnableIf
{
} ;

template < typename ReturnType_ >
struct EnableIf< true, ReturnType_ >
{
	typedef ReturnType_ ReturnType ;
} ;

template < bool Condition, typename ReturnType_ = void >
struct DisableIf
{
} ;

template < typename ReturnType_ >
struct DisableIf< false, ReturnType_ >
{
	typedef ReturnType_ ReturnType ;
} ;

// Warning : HasXXX< T > does NOT work for reference typedefs !
#define BOGUS_DEFINE_HAS_TYPE( TypeName ) \
	template < typename BaseType > 	      \
	struct Has##TypeName	              \
	{	                                  \
	private:                              \
		enum { True = 1, False = 2 } ;    \
		typedef char  TrueType[  True ] ; \
		typedef char FalseType[ False ] ; \
										  \
		template< typename T >			  \
		static const  TrueType& check( const typename T::TypeName* ) ; \
		template< typename >			  \
		static const FalseType& check( ... ) ; \
	public:								  \
		enum { Value = ( True == sizeof( check< BaseType >( 0 ) ) ) } ;\
	}

BOGUS_DEFINE_HAS_TYPE( ReturnType ) ;
BOGUS_DEFINE_HAS_TYPE( ConstTransposeReturnType ) ;
BOGUS_DEFINE_HAS_TYPE( Base ) ;

// Static assertions

template < bool Assertion >
struct StaticAssert
{
	enum {
		BLOCKS_MUST_BE_SQUARE_OR_HAVE_DYNAMIC_DIMENSIONS,
		BLOCKS_MUST_HAVE_FIXED_DIMENSIONS,
		MATRICES_ORDERING_IS_INCONSISTENT,
		TRANSPOSE_MAKES_NO_SENSE_IN_THIS_CONTEXT,
		TRANSPOSE_IS_NOT_DEFINED_FOR_THIS_BLOCK_TYPE,
		OPERANDS_HAVE_INCONSISTENT_FLAGS,
		UNORDERED_INSERTION_WITH_COMPRESSED_INDEX
	} ;
} ;

template < >
struct StaticAssert< false >
{
} ;

#define BOGUS_STATIC_ASSERT( test, message ) (void) ::bogus::StaticAssert< test >::message

//! Const mapped array, used for Mapped Block Matrices
template< typename Element >
class ConstMappedArray
{
public:

	typedef ConstMappedArray< Element > Type ;
	enum { is_mutable = 0 } ;

	ConstMappedArray ( )
		: m_data( 0 ), m_size( 0 )
	{}

	ConstMappedArray ( const Element* data, std::size_t size )
		: m_data( data ), m_size( size )
	{}

	void setData( const Element* data, std::size_t size )
	{
		m_data = data ;
		m_size = size ;
	}

	const Element* data() const { return m_data ; }

	inline std::size_t size() const { return m_size ; }
	inline bool empty() const { return 0 == m_size ; }

	const Element* begin() const { return data() ; }
	const Element* end() const { return data() + size() ; }

	const Element& operator[]( std::size_t idx ) const
	{ return m_data[idx] ; }

	inline void resize( std::size_t s)
	{ assert( !m_data || m_size == s ) ; (void) s ; }
	inline void reserve( std::size_t )
	{ }
	inline void clear( )
	{ }
	inline void assign( std::size_t s, const Element&)
	{ assert( !m_data || m_size == s ) ; (void) s ; }

private:
	const Element* m_data ;
	std::size_t    m_size ;
} ;

} //namespace bogus


#endif



