/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#ifndef BOGUS_BLOCK_SOLVERS_FWD_HPP
#define BOGUS_BLOCK_SOLVERS_FWD_HPP

namespace bogus {

//! Configuration properties for Krylov solvers
namespace krylov
{
enum Method
{
	CG,			//!< Conjugate Gradient. \sa krylov::solvers::CG
	BiCG,			//!< BiConjugate Gradient \sa krylov::solvers::BiCG
	BiCGSTAB, 		//!< BiConjugate Gradient Stabilized \sa krylov::solvers::BiCGSTAB
	CGS, 			//!< Conjugate Gradient Squared \sa krylov::solvers::CGS
	GMRES,			//!< Generalized Minimal Residual \sa krylov::solvers::GMRES
	TFQMR			//!< Tranpose-free Quasi Minimal Residual \sa krylov::solvers::TFQMR
} ;

} // namespace iterative_linear_solvers

template < typename MatrixType >
struct ProblemTraits ;

template < template <typename> class Method, typename BlockMatrixType >
class ConstrainedSolverBase ;

template < typename BlockMatrixType >
class GaussSeidel ;

template < typename BlockMatrixType >
class ProjectedGradient ;

template < typename MatrixType >
class TrivialPreconditioner ;

template < typename BlockMatrixType,
		   template< typename BlockMatrixT > class PreconditionerType = TrivialPreconditioner >
class Krylov ;

}

#endif
