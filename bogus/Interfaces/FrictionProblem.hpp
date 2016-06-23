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

#ifndef BOGUS_FRICTION_PROBLEM_HPP
#define BOGUS_FRICTION_PROBLEM_HPP

#include "../Core/Block.hpp"
#include "../Core/BlockSolvers.fwd.hpp"

#include "../Extra/SecondOrder.fwd.hpp"

#include "../Core/Utils/Signal.hpp"
#include "../Core/ExternalForce.hh"

#include "../../StrandSim/Core/Definitions.hh"


namespace bogus
{

template< unsigned Dimension >
struct PrimalFrictionProblem
{
	// Primal Data
	//! M -- mass matrix
	SparseBlockMatrix< Eigen::MatrixXd > M ;
	//! E -- local rotation matrix ( world <-> contact basis )
	bogus::SparseBlockMatrix< Eigen::Matrix< double, Dimension, Dimension > > E ;

	typedef Eigen::Matrix< double, Dimension, Eigen::Dynamic > HBlock ;
	//! H -- deformation gradient ( generalized coordinates <-> 3D world )
	SparseBlockMatrix< HBlock, UNCOMPRESSED > H;

	//! External forces
	 Eigen::VectorXd f ;
	//! Free velocity ( such that u = Hv + w )
	 Eigen::VectorXd w ;
	//! Coulomb friction coefficients
	 Eigen::VectorXd mu ;

	// Cached data

	//! M^-1
	SparseBlockMatrix< LU< Eigen::MatrixBase< Eigen::MatrixXd > > > MInv ;


    /*! Returns true if at least one force has been recomputed
      ( and therefore the velocities are note up to date )
    */
    bool updateExternalForces( Eigen::VectorXd& m_v, Eigen::VectorXd m_r, const std::vector < unsigned >& m_startDofs );

    void updateObjectLHS ( unsigned oId, Eigen::MatrixXd &M, std::vector<unsigned> ndof ) ;

    bool hasExternalForces() const { return ! m_externalForces.empty() ; }
    std::vector< ExternalForce * > m_externalForces ;

} ;


template< unsigned Dimension >
struct DualFrictionProblem
{
	typedef SparseBlockMatrix< Eigen::Matrix< double, Dimension, Dimension, Eigen::RowMajor >,
							   SYMMETRIC > WType ;

	typedef GaussSeidel< WType > GaussSeidelType ;
	typedef ProjectedGradient< WType > ProjectedGradientType ;

	typedef SOCLaw< Dimension, double, true  > CoulombLawType	;
	typedef SOCLaw< Dimension, double, false > SOCLawType	;

	//! W -- Delassus operator
	WType W ;

	//! Rhs ( such that u = Wr + b )
	Eigen::VectorXd b ;

	//! Coulomb friction coefficients
	Eigen::VectorXd mu ;

	//! Computes this DualFrictionProblem from the given \p primal
	void computeFrom( PrimalFrictionProblem< Dimension >& primal ) ;

	//! Solves this problem
	/*!
	  \param gs The GaussSeidel< WType > solver to use
	  \param r  Both the initial guess and the result
	  \param staticProblem If true, solve this problem as a \b SOCQP instead of a Coulomb Friction problem
	  \returns the error as returned by the GaussSeidel::solve() function
	  */
	double solveWith( GaussSeidelType &gs, double * r, const bool staticProblem = false ) const ;
	double solveWith( ProjectedGradientType &pg, double * r ) const ;

	//! Evaluate a residual using the GS's error function
	/*!
	  \param gs The GaussSeidel< WType > solver to use
	  \param r  Both the current force
	  \param staticProblem If true, eval this problem as a \b SOCQP instead of a Coulomb Friction problem

	  \returns the error as returned by the GaussSeidel::eval() function
	  */
	double evalWith( const GaussSeidelType &gs, const double * r, const bool staticProblem = false ) const ;
	double evalWith( const ProjectedGradientType &gs, const double * r ) const ;

	//! Solves this problem using the Cadoux algorithm ( with fixed-point iteration )
	/*!
	  See \cite ACLM11
	  \param gs The GaussSeidel< WType > solver to use
	  \param r  Both the initial guess and the result
	  \param fpIterations Number of fixed-point iterations
	  \param callback 0, or a pointer to a user-defined function that takes ( unsigned iteration, double residual ) as arguments
	  \returns the error as returned by the GaussSeidel::solve() function
	  */
	double solveCadoux( GaussSeidelType &gs, double * r, const unsigned fpIterations,
		   const Signal< unsigned, double >* callback = 0 ) const ;
	double solveCadoux( ProjectedGradientType &pg, double * r, const unsigned fpIterations,
		   const Signal< unsigned, double >* callback = 0 ) const ;


	//! \warning To use the permutation releated functions, all the blocks have to have the same size
	void applyPermutation( const std::vector< std::size_t > &permutation ) ;
	void undoPermutation() ;
	bool permuted() const { return !m_permutation.empty() ; }

	const std::vector< std::size_t > &permutation() const { return m_permutation ; }
	const std::vector< std::size_t > &invPermutation() const { return m_invPermutation ; }
private:

	// Current permutation of contact indices
	std::vector< std::size_t > m_permutation ;
	std::vector< std::size_t > m_invPermutation ;
} ;

} //namespace bogus

#endif
