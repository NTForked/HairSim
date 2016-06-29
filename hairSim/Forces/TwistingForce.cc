#include "TwistingForce.hh"
#include "ViscousOrNotViscous.hh"
#include "../Math/BandMatrix.h"

template<typename ViscousT>
Scalar TwistingForce<ViscousT>::localEnergy( const ElasticStrand& strand, StrandState& geometry,
        const IndexType vtx )
{
    const Scalar kt = ViscousT::kt( strand, vtx );
    const Scalar undefTwist = ViscousT::thetaBar( strand, vtx );
    const Scalar ilen = strand.m_invVoronoiLengths[vtx];
    const Scalar twist = geometry.m_twists[vtx];

    return 0.5 * kt * square( twist - undefTwist ) * ilen;
}

template<typename ViscousT>
void TwistingForce<ViscousT>::computeLocal( TwistingForce::LocalForceType& localF,
        const ElasticStrand& strand, StrandState& geometry, const IndexType vtx )
{
    const Scalar kt = ViscousT::kt( strand, vtx );
    const Scalar undefTwist = ViscousT::thetaBar( strand, vtx );
    const Scalar ilen = strand.m_invVoronoiLengths[vtx];
    const Scalar twist = geometry.m_twists[vtx];

    localF = -kt * ilen * ( twist - undefTwist ) * geometry.m_gradTwists[vtx];
}

template<typename ViscousT>
void TwistingForce<ViscousT>::computeLocal( TwistingForce::LocalThetaForceType& localF,
        const ElasticStrand& strand, StrandState& geometry, const IndexType vtx )
{
    const Scalar kt = ViscousT::kt( strand, vtx );
    const Scalar undefTwist = ViscousT::thetaBar( strand, vtx );
    const Scalar ilen = strand.m_invVoronoiLengths[vtx];
    const Scalar twist = geometry.m_twists[vtx];
    const Vec11x& gradTwist = geometry.m_gradTwists[vtx];
    Vec2 thetaGradTwist( gradTwist[3], gradTwist[7] );

    localF = -kt * ilen * ( twist - undefTwist ) * thetaGradTwist;
}

template<typename ViscousT>
void TwistingForce<ViscousT>::computeLocal( TwistingForce::LocalJacobianType& localJ,
        const ElasticStrand& strand, StrandState& geometry, const IndexType vtx )
{
    const Scalar kt = ViscousT::kt( strand, vtx );
    const Scalar ilen = strand.m_invVoronoiLengths[vtx];
    const Mat11x& gradTwistSquared = geometry.m_gradTwistsSquared[vtx];

    localJ = -kt * ilen * gradTwistSquared ;
    if( strand.requiresExactJacobian() )
    {
        const Scalar undeformedTwist = ViscousT::thetaBar( strand, vtx );
        const Scalar twist = geometry.m_twists[vtx];
        const Mat11x& hessTwist = geometry.m_hessTwists[vtx];
        localJ += -kt * ilen * ( twist - undeformedTwist ) * hessTwist ;
    }
}

template<typename ViscousT>
void TwistingForce<ViscousT>::computeLocal( TwistingForce::LocalThetaJacobianType& localJ,
        const ElasticStrand& strand, StrandState& geometry, const IndexType vtx )
{
    const Scalar kt = ViscousT::kt( strand, vtx );
    const Scalar undeformedTwist = ViscousT::thetaBar( strand, vtx );
    const Scalar ilen = strand.m_invVoronoiLengths[vtx];
    const Vec11x& gradTwist = geometry.m_gradTwists[vtx];
    Vec2 thetaGradTwist( gradTwist[3], gradTwist[7] );

    // There is no twist Hessian on the theta coordinate.
    localJ = -kt * ilen * ( thetaGradTwist * thetaGradTwist.transpose() );
}

template<typename ViscousT>
void TwistingForce<ViscousT>::addInPosition( VecXx& globalForce, const IndexType vtx,
        const LocalForceType& localForce )
{
    globalForce.segment<11>( 4 * ( vtx - 1 ) ) += localForce;
}

template<typename ViscousT>
void TwistingForce<ViscousT>::addInPosition( VecXx& globalForce, const IndexType vtx,
        const LocalThetaForceType& localForce )
{
    globalForce.segment<2>( vtx - 1 ) += localForce;
}

template<typename ViscousT>
void TwistingForce<ViscousT>::addInPosition( JacobianMatrixType& globalJacobian,
        const IndexType vtx, const LocalJacobianType& localJacobian )
{
    globalJacobian.localStencilAdd<11>( 4 * ( vtx - 1 ), localJacobian );
}

template<typename ViscousT>
void TwistingForce<ViscousT>::addInPosition( TriDiagonalMatrixType& globalJacobian,
        const IndexType vtx, const LocalThetaJacobianType& localJacobian )
{
    globalJacobian.localStencilAdd<2>( vtx - 1, localJacobian );
}

template<typename ViscousT>
void TwistingForce<ViscousT>::accumulateCurrentE( Scalar& energy, ElasticStrand& strand )
{
    for ( IndexType vtx = s_first; vtx < strand.getNumVertices() - s_last; ++vtx )
    {
        energy += localEnergy( strand, strand.getCurrentState(), vtx );
    }
}

template<typename ViscousT>
void TwistingForce<ViscousT>::accumulateCurrentF( VecXx& force, ElasticStrand& strand )
{
    for ( IndexType vtx = s_first; vtx < strand.getNumVertices() - s_last; ++vtx )
    {
        LocalForceType localF;
        computeLocal( localF, strand, strand.getCurrentState(), vtx );
        addInPosition( force, vtx, localF );
    }
}

template<typename ViscousT>
void TwistingForce<ViscousT>::accumulateCurrentJ( JacobianMatrixType& Jacobian,
        ElasticStrand& strand )
{
    for ( IndexType vtx = s_first; vtx < strand.getNumVertices() - s_last; ++vtx )
    {
        LocalJacobianType localJ;
        computeLocal( localJ, strand, strand.getCurrentState(), vtx );
        addInPosition( Jacobian, vtx, localJ );
    }
}

template class TwistingForce<NonViscous> ;
template class TwistingForce<Viscous> ;
