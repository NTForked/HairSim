#include "BendingForce.hh"
#include "ViscousOrNotViscous.hh"
#include "../Core/ElasticStrandUtils.hh"
#include "../Math/BandMatrix.h"

template<typename ViscousT>
Scalar BendingForce<ViscousT>::localEnergy( const ElasticStrand& strand, StrandState& geometry,
        const IndexType vtx )
{
    const Mat2x& B = ViscousT::bendingMatrix( strand, vtx );
    const Vec2& kappaBar = ViscousT::kappaBar( strand, vtx );
    const Scalar ilen = strand.m_invVoronoiLengths[vtx];
    const Vec2& kappa = geometry.m_kappas[vtx];

    return 0.5 * ilen * ( kappa - kappaBar ).dot( Vec2( B * ( kappa - kappaBar ) ) );
}

template<typename ViscousT>
void BendingForce<ViscousT>::computeLocal( BendingForce::LocalForceType& localF,
        const ElasticStrand& strand, StrandState& geometry, const IndexType vtx )
{
    const Mat2x& B = ViscousT::bendingMatrix( strand, vtx );
    const Vec2& kappaBar = ViscousT::kappaBar( strand, vtx );
    const Scalar ilen = strand.m_invVoronoiLengths[vtx];
    const Vec2& kappa = geometry.m_kappas[vtx];
    const GradKType& gradKappa = geometry.m_gradKappas[vtx];

    localF = -ilen * gradKappa * B * ( kappa - kappaBar );
}

template<typename ViscousT>
void BendingForce<ViscousT>::computeLocal( BendingForce::LocalThetaForceType& localF,
        const ElasticStrand& strand, StrandState& geometry, const IndexType vtx )
{
    const GradKType& gradKappa = geometry.m_gradKappas[vtx];
    Mat2x thetaGradKappa;
    thetaGradKappa( 0, 0 ) = gradKappa( 3, 0 );
    thetaGradKappa( 0, 1 ) = gradKappa( 3, 1 );
    thetaGradKappa( 1, 0 ) = gradKappa( 7, 0 );
    thetaGradKappa( 1, 1 ) = gradKappa( 7, 1 );

    const Mat2x& B = ViscousT::bendingMatrix( strand, vtx );
    const Vec2& kappaBar = ViscousT::kappaBar( strand, vtx );
    const Scalar ilen = strand.m_invVoronoiLengths[vtx];
    const Vec2& kappa = geometry.m_kappas[vtx];

    localF = -ilen * thetaGradKappa * B * ( kappa - kappaBar );
}

template<typename ViscousT>
void BendingForce<ViscousT>::computeLocal( BendingForce::LocalJacobianType& localJ,
        const ElasticStrand& strand, StrandState& geometry, const IndexType vtx )
{
    localJ = geometry.m_bendingProducts[vtx];

    if ( strand.requiresExactJacobian() )
    {
        const Mat2x& bendingMatrixBase = strand.m_parameters.bendingMatrixBase();
        const Vec2& kappaBar = ViscousT::kappaBar( strand, vtx );
        const Vec2& kappa = geometry.m_kappas[vtx];
        const std::pair<LocalJacobianType, LocalJacobianType>& hessKappa =
                geometry.m_hessKappas[vtx];
        const Vec2& temp = bendingMatrixBase * ( kappa - kappaBar );

        localJ += temp( 0 ) * hessKappa.first + temp( 1 ) * hessKappa.second;
    }

    const Scalar ilen = strand.m_invVoronoiLengths[vtx];
    localJ *= -ilen * ViscousT::bendingCoefficient( strand, vtx );
}

template<typename ViscousT>
void BendingForce<ViscousT>::computeLocal( BendingForce::LocalThetaJacobianType& localJ,
        const ElasticStrand& strand, StrandState& geometry, const IndexType vtx )
{
    const GradKType& gradKappa = geometry.m_gradKappas[vtx];
    Mat2x thetaGradKappa;
    thetaGradKappa( 0, 0 ) = gradKappa( 3, 0 );
    thetaGradKappa( 0, 1 ) = gradKappa( 3, 1 );
    thetaGradKappa( 1, 0 ) = gradKappa( 7, 0 );
    thetaGradKappa( 1, 1 ) = gradKappa( 7, 1 );

    const Mat2x& bendingMatrixBase = strand.m_parameters.bendingMatrixBase();
    symBProduct<2>( localJ, bendingMatrixBase, thetaGradKappa );

    if ( strand.requiresExactJacobian() )
    {
        const Vec2& kappaBar = ViscousT::kappaBar( strand, vtx );
        const Vec2& kappa = geometry.m_kappas[vtx];
        const std::pair<LocalThetaJacobianType, LocalThetaJacobianType>& hessKappa =
                geometry.m_thetaHessKappas[vtx];
        const Vec2& temp = bendingMatrixBase * ( kappa - kappaBar );
        localJ += temp( 0 ) * hessKappa.first + temp( 1 ) * hessKappa.second;
    }

    const Scalar ilen = strand.m_invVoronoiLengths[vtx];
    localJ *= -ilen * ViscousT::bendingCoefficient( strand, vtx );
}

template<typename ViscousT>
void BendingForce<ViscousT>::addInPosition( VecXx& globalForce, const IndexType vtx,
        const LocalForceType& localForce )
{
    globalForce.segment<11>( 4 * ( vtx - 1 ) ) += localForce;
}

template<typename ViscousT>
void BendingForce<ViscousT>::addInPosition( VecXx& globalForce, const IndexType vtx,
        const LocalThetaForceType& localForce )
{
    globalForce.segment<2>( vtx - 1 ) += localForce;
}

template<typename ViscousT>
void BendingForce<ViscousT>::addInPosition( JacobianMatrixType& globalJacobian, const IndexType vtx,
        const LocalJacobianType& localJacobian )
{
    globalJacobian.localStencilAdd<11>( 4 * ( vtx - 1 ), localJacobian );
}

template<typename ViscousT>
void BendingForce<ViscousT>::addInPosition( TriDiagonalMatrixType& globalJacobian,
        const IndexType vtx, const LocalThetaJacobianType& localJacobian )
{
    globalJacobian.localStencilAdd<2>( vtx - 1, localJacobian );
}

template<typename ViscousT>
void BendingForce<ViscousT>::accumulateCurrentE( Scalar& energy, ElasticStrand& strand )
{
    for ( IndexType vtx = s_first; vtx < strand.getNumVertices() - s_last; ++vtx )
    {
        energy += localEnergy( strand, strand.getCurrentState(), vtx );
    }
}

template<typename ViscousT>
void BendingForce<ViscousT>::accumulateCurrentF( VecXx& force, ElasticStrand& strand )
{
    for ( IndexType vtx = s_first; vtx < strand.getNumVertices() - s_last; ++vtx )
    {
        LocalForceType localF;
        computeLocal( localF, strand, strand.getCurrentState(), vtx );
        addInPosition( force, vtx, localF );
    }
}

template<typename ViscousT>
void BendingForce<ViscousT>::accumulateCurrentJ( JacobianMatrixType& Jacobian,
        ElasticStrand& strand )
{
    for ( IndexType vtx = s_first; vtx < strand.getNumVertices() - s_last; ++vtx )
    {
        LocalJacobianType localJ;
        computeLocal( localJ, strand, strand.getCurrentState(), vtx );
        addInPosition( Jacobian, vtx, localJ );
    }
}

template class BendingForce<NonViscous>;
template class BendingForce<Viscous>;
