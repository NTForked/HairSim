/*
 * TwistingForce.hh
 *
 *  Created on: 12/07/2011
 *      Author: Jean-Marie Aubry <jaubry@wetafx.co.nz>
 */

#ifndef TWISTINGFORCE_HH_
#define TWISTINGFORCE_HH_

#include "ViscousOrNotViscous.hh"
#include "../Core/ElasticStrand.hh"

namespace strandsim
{

template<typename ViscousT = NonViscous>
class TwistingForce: public ForceBase
{
public:
    TwistingForce()
    {}
    virtual ~TwistingForce()
    {}

public:
    static const IndexType s_first = 1; // The first index on which this force can apply
    static const IndexType s_last = 1; // The last index (counting from the end)

    typedef Eigen::Matrix<Scalar, 11, 1> LocalForceType;
    typedef Vec2x LocalThetaForceType;
    typedef Eigen::Matrix<Scalar, 11, 11> LocalJacobianType;
    typedef Mat2x LocalThetaJacobianType;

    static std::string getName()
    {
        return ViscousT::getName() + "twisting";
    }

    static Scalar localEnergy( const ElasticStrand& strand, StrandState& geometry,
            const IndexType vtx );

    static void computeLocal( LocalForceType& localF, const ElasticStrand& strand,
            StrandState& geometry, const IndexType vtx );
    static void computeLocal( LocalThetaForceType& localF, const ElasticStrand& strand,
            StrandState& geometry, const IndexType vtx );

    static void computeLocal( LocalJacobianType& localJ, const ElasticStrand& strand,
            StrandState& geometry, const IndexType vtx );
    static void computeLocal( LocalThetaJacobianType& localJ, const ElasticStrand& strand,
            StrandState& geometry, const IndexType vtx );

    static void addInPosition( VecXx& globalForce, const IndexType vtx,
            const LocalForceType& localForce );
    static void addInPosition( VecXx& globalForce, const IndexType vtx,
            const LocalThetaForceType& localForce );

    static void addInPosition( JacobianMatrixType& globalJacobian, const IndexType vtx,
            const LocalJacobianType& localJacobian );
    static void addInPosition( TriDiagonalMatrixType& globalJacobian, const IndexType vtx,
            const LocalThetaJacobianType& localJacobian );

    static void accumulateCurrentE( Scalar& energy, ElasticStrand& strand );
    static void accumulateCurrentF( VecXx& force, ElasticStrand& strand );
    static void accumulateCurrentJ( JacobianMatrixType& Jacobian, ElasticStrand& strand );
};

}

#endif /* TWISTINGFORCE_HH_ */
