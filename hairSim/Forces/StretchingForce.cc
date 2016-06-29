#include "StretchingForce.hh"
#include "ViscousOrNotViscous.hh"
#include "../Math/BandMatrix.h"

template<typename ViscousT>
Scalar StretchingForce<ViscousT>::localEnergy( const ElasticStrand& strand, StrandState& geometry,
        const IndexType vtx )
{
    const Scalar ks = ViscousT::ks( strand, vtx );
    const Scalar restLength = ViscousT::ellBar( strand, vtx );

    const Scalar length = geometry.m_lengths[vtx];

    return 0.5 * ks * square( length / restLength - 1.0 ) * restLength;
}

template<typename ViscousT>
void StretchingForce<ViscousT>::computeLocal( StretchingForce::LocalForceType& localF,
        const ElasticStrand& strand, StrandState& geometry, const IndexType vtx )
{
    const Scalar ks = ViscousT::ks( strand, vtx );
    const Scalar restLength = ViscousT::ellBar( strand, vtx );

    const Scalar length = geometry.m_lengths[vtx];
    const Vec3& edge = geometry.m_tangents[ vtx ];

    Vec3 f = ks * ( length / restLength - 1.0 ) * edge ;

    // if( strand.getGlobalIndex() == 1 ) std::cout << "{"<< strand.getGlobalIndex() <<"}{"<< vtx <<"} l: " << length << " rl: " << restLength << " f: " << f.norm() << std::endl;

    localF.segment<3>( 0 ) = f;
    localF.segment<3>( 3 ) = -f;
}

template<typename ViscousT>
void StretchingForce<ViscousT>::computeLocal( StretchingForce::LocalJacobianType& localJ,
        const ElasticStrand& strand, StrandState& geometry, const IndexType vtx )
{
    const Scalar ks = ViscousT::ks( strand, vtx );
    const Scalar restLength = ViscousT::ellBar( strand, vtx );

    const Scalar length = geometry.m_lengths[vtx];
    const Vec3& edge = geometry.m_tangents[ vtx ];

    bool useApprox = !strand.requiresExactJacobian() && length < restLength ;

    Mat3x M ;
    if( useApprox )
    {
        M = ks / restLength * ( edge * edge.transpose() ) ;
    } else {
        M = ks
                * ( ( 1.0 / restLength - 1.0 / length ) * Mat3x::Identity()
                    + 1.0 / length * ( edge * edge.transpose() ) );
    }

    localJ.block<3, 3>( 0, 0 ) = localJ.block<3, 3>( 3, 3 ) = -M;
    localJ.block<3, 3>( 0, 3 ) = localJ.block<3, 3>( 3, 0 ) = M;
}

template<typename ViscousT>
void StretchingForce<ViscousT>::addInPosition( ForceVectorType& globalForce, const IndexType vtx,
        const LocalForceType& localForce )
{
    globalForce.segment<3>( 4 * vtx ) += localForce.segment<3>( 0 );
    globalForce.segment<3>( 4 * ( vtx + 1 ) ) += localForce.segment<3>( 3 );
}

template<typename ViscousT>
void StretchingForce<ViscousT>::addInPosition( JacobianMatrixType& globalJacobian,
        const IndexType vtx, const LocalJacobianType& localJacobian )
{
    globalJacobian.edgeStencilAdd<6>( 4 * vtx, localJacobian );
}

template<typename ViscousT>
void StretchingForce<ViscousT>::accumulateCurrentE( Scalar& energy, ElasticStrand& strand )
{
    for ( IndexType vtx = s_first; vtx < strand.getNumVertices() - s_last; ++vtx )
    {
        energy += localEnergy( strand, strand.getCurrentState(), vtx );
    }
}

template<typename ViscousT>
void StretchingForce<ViscousT>::accumulateCurrentF( VecXx& force, ElasticStrand& strand )
{
    for ( IndexType vtx = s_first; vtx < strand.getNumVertices() - s_last; ++vtx )
    {
        LocalForceType localF;
        computeLocal( localF, strand, strand.getCurrentState(), vtx );
        addInPosition( force, vtx, localF );
    }
}

template<typename ViscousT>
void StretchingForce<ViscousT>::accumulateCurrentJ( JacobianMatrixType& Jacobian,
        ElasticStrand& strand )
{
    for ( IndexType vtx = s_first; vtx < strand.getNumVertices() - s_last; ++vtx )
    {
        LocalJacobianType localJ;
        computeLocal( localJ, strand, strand.getCurrentState(), vtx );
        addInPosition( Jacobian, vtx, localJ );
    }
}


template class StretchingForce<NonViscous> ;
template class StretchingForce<Viscous> ;
