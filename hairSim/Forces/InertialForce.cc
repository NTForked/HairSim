#include "InertialForce.hh"
#include "../Core/ElasticStrand.hh"
#include "../Dynamic/StrandDynamicTraits.hh"
#include "../Core/ElasticStrandUtils.hh"
#include "../Core/BandMatrix.hh"

namespace strandsim {

Mat3x InertialForce::s_Omega_Cross = Mat3x::Zero() ;
Mat3x InertialForce::s_Omega_Cross_Cross = Mat3x::Zero() ;
Mat3x InertialForce::s_dOmega_dt_Cross = Mat3x::Zero() ;
Vec3x InertialForce::s_accOrigin = Vec3x::Zero() ;

InertialForce::InertialForce()
{
}

InertialForce::~InertialForce()
{
}

void InertialForce::computeLocal( LocalForceType& localF, const ElasticStrand& strand,
        const StrandState& geometry, const IndexType vtx )
{
    const Scalar m = strand.getVertexMass( vtx );
    const Vec3x& P = geometry.getVertex( vtx ) ;
    const Vec3x vr = strand.dynamics().getDisplacement( vtx ) / strand.getParameters().getDt() ;

    const Vec3x& ae = s_accOrigin + s_dOmega_dt_Cross * P
            + s_Omega_Cross_Cross * P ;
            + 2. * s_Omega_Cross * vr ;

    localF = - m * ae ;
}

void InertialForce::computeLocal( LocalJacobianType& localJ, const ElasticStrand& strand,
        const StrandState& geometry, const IndexType vtx )
{
    const Scalar m = strand.getVertexMass( vtx );

    localJ = - m * ( s_dOmega_dt_Cross + s_Omega_Cross_Cross
                    + 2. * s_Omega_Cross / strand.getParameters().getDt()
                     ) ;
}


void InertialForce::addInPosition( ForceVectorType& globalForce, const IndexType vtx,
        const LocalForceType& localForce )
{
    globalForce.segment<3>( 4 * vtx ) += localForce;
}

void InertialForce::addInPosition( JacobianMatrixType& globalJacobian, const IndexType vtx,
        const LocalJacobianType& localJacobian )
{
    globalJacobian.localStencilAdd<3>( 4 * vtx, localJacobian );
}

void InertialForce::setFrameAccelerations(
            const Vec3x& Omega, const Vec3x& dOmega_dt,
            const Vec3x& accOrigin )
{
    s_Omega_Cross = crossMat( Omega ) ;
    s_Omega_Cross_Cross = s_Omega_Cross * s_Omega_Cross ;
    s_dOmega_dt_Cross = crossMat( dOmega_dt ) ;
    s_accOrigin = accOrigin ;
}

} // namespace strandsim
