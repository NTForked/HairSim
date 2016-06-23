/*
 * GravitationForce.cc
 *
 *  Created on: 14/07/2011
 *      Author: Jean-Marie Aubry <jaubry@wetafx.co.nz>
 */

#include "GravitationForce.hh"
#include "../Core/ElasticStrand.hh"
#include "../Core/BandMatrix.hh"

GravitationForce::GravitationForce()
{
    // TODO Auto-generated constructor stub
}

GravitationForce::~GravitationForce()
{
    // TODO Auto-generated destructor stub
}

Scalar GravitationForce::localEnergy( const ElasticStrand& strand, const StrandState& geometry,
        const IndexType vtx )
{
    return -strand.m_vertexMasses[vtx] * geometry.getVertex( vtx ).dot( s_gravity );
}

template<>
void GravitationForce::computeLocal<Scalar>( Scalar& localE, const ElasticStrand& strand,
        const StrandState& geometry, const IndexType vtx )
{
    localE = -strand.m_vertexMasses[vtx] * geometry.getVertex( vtx ).dot( s_gravity );
}

template<>
void GravitationForce::computeLocal<Vec3x>( Vec3x& localF, const ElasticStrand& strand,
        const StrandState& geometry, const IndexType vtx )
{
    localF = strand.m_vertexMasses[vtx] * s_gravity;
}

template<>
void GravitationForce::computeLocal<Mat3x>( Mat3x& localJ, const ElasticStrand& strand,
        const StrandState& geometry, const IndexType vtx )
{
    localJ.setZero(); // Jacobian is zero
}

template<>
void GravitationForce::addInPosition<VecXx, Vec3x>( VecXx& globalForce, const IndexType vtx,
        const Vec3x& localForce )
{
    globalForce.segment<3>( 4 * vtx ) += localForce;
}

template<>
void GravitationForce::addInPosition<JacobianMatrixType, Mat3x>( JacobianMatrixType& globalJacobian,
        const IndexType vtx, const Mat3x& localJacobian )
{
    globalJacobian.localStencilAdd<3>( 4 * vtx, localJacobian );
}

Vec3x GravitationForce::s_gravity( 0.0, -981.0, 0.0 ); // Acceleration vector, in cm/s^2 to match Maya's units. Also, the y-axis is vertical.
