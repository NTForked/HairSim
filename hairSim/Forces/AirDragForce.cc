#include "AirDragForce.hh"
#include "../Core/ElasticStrand.hh"
#include "../Core/ElasticStrandUtils.hh"
#include "../Dynamic/StrandDynamicTraits.hh"
#include "../Core/BandMatrix.hh"

Vec3x AirDragForce::s_velOrigin = Vec3x::Zero() ;
Mat3x AirDragForce::s_Omega_Cross = Mat3x::Zero() ;

AirDragForce::AirDragForce()
{}

AirDragForce::~AirDragForce()
{}

void AirDragForce::computeLocal( LocalForceType& localF, const ElasticStrand& strand,
        const StrandState& geometry, const IndexType vtx )
{
    const Scalar nu = strand.getParameters().getAirDrag() ;
    const Scalar len = strand.m_VoronoiLengths[ vtx ];

    const Vec3x ve = s_velOrigin + s_Omega_Cross * geometry.getVertex( vtx ) ;
    const Vec3x vr = strand.dynamics().getDisplacement( vtx ) / strand.getParameters().getDt() ;

    localF = - len * nu *  ( ve + vr ) ;
}

void AirDragForce::computeLocal( LocalJacobianType& localJ, const ElasticStrand& strand,
        const StrandState& geometry, const IndexType vtx )
{
    const Scalar nu = strand.getParameters().getAirDrag() ;
    const Scalar len = strand.m_VoronoiLengths[ vtx ];

    localJ = -len * nu * ( Mat3x::Identity() / strand.getParameters().getDt() + s_Omega_Cross );
}

void AirDragForce::addInPosition( ForceVectorType& globalForce, const IndexType vtx,
        const LocalForceType& localForce )
{
    globalForce.segment<3>( 4 * vtx ) += localForce;
}

void AirDragForce::addInPosition( JacobianMatrixType& globalJacobian, const IndexType vtx,
        const LocalJacobianType& localJacobian )
{
    globalJacobian.localStencilAdd<3>( 4 * vtx, localJacobian );
}

void AirDragForce::setFrameVelocities(
            const Vec3x& Omega,
            const Vec3x& velOrigin )
{
    s_velOrigin = velOrigin ;
    s_Omega_Cross = crossMat( Omega ) ;
}
