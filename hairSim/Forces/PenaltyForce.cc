#include "PenaltyForce.hh"

#include "../Core/ElasticStrand.hh"
#include "../Core/BandMatrix.hh"

PenaltyForce::PenaltyForce()
    : TunneledBandRuntimeForce(),
      m_stiffness( 0.0 ),
      m_thickness( 0.0 )
{}

PenaltyForce::~PenaltyForce()
{}

void PenaltyForce::setParameters( const std::vector< ElasticStrand* > &strands, const Vec3x& origin, const Scalar stiffness,
                                 const Scalar restLength, const bool allowCompression )
{
    m_origin = origin ;
    m_stiffness = stiffness ;
    m_restLength = restLength ;
    m_allowCompression = allowCompression ;

    for( auto it = m_vertices.begin() ; it != m_vertices.end() ; it += VerticesSet::nextSet )
    {
        if( (*it).first < (int) strands.size() && strands[ (*it).first ] )
        {
            strands[ (*it).first ]->invalidatePhysics() ;
        }
    }
}

bool PenaltyForce::breakSprings( const std::vector< ElasticStrand* > &strands, const Scalar breakingForce )
{
    Vec3x localF ;
    for( auto vIt = m_vertices.begin() ; vIt != m_vertices.end() ; )
    {
       localF.setZero() ;
       ElasticStrand& strand = *strands[ (*vIt).first ] ;
       const Vec3x &v = strand.getVertex( (*vIt).second ) ;

       computeLocalF( v, strand, localF ) ;
       if( localF.norm() > breakingForce )
       {
           // Clean up m_sets if necessary
           if( vIt.getSet().size() == 1 )
           {
               strand.removeRuntimeForce( this ) ;
               m_sets.erase( &strand ) ;
           }

           vIt = m_vertices.erase( vIt ) ;

       } else {
           ++ vIt ;
       }
    }

    return ! m_vertices.empty() ;

}

void PenaltyForce::accumulateEFJ( StrandState& geometry, const ElasticStrand& strand )
{
    auto setIt = m_sets.find( &strand ) ;
    if ( setIt == m_sets.end() ) return ;

    const VerticesSet::IndexedSetType::SetType &strandVertices = setIt->second.getSet() ;

    Vec3x localF ;
    Mat3x localJ ;
    for( auto vIt = strandVertices.begin() ; vIt != strandVertices.end() ; ++vIt )
    {
       const Vec3x &v = geometry.getVertex( *vIt ) ;
       geometry.m_totalEnergy += computeLocalE( v, strand ) ;
       computeLocalF( v, strand, localF ) ;
       geometry.m_totalForce.segment< 3 >( 4 * *vIt ) += localF ;
       computeLocalJ( v, strand, localJ ) ;
       geometry.m_totalJacobian->localStencilAdd< 3 >( 4 * *vIt, localJ ) ;
    }

}

void PenaltyForce::accumulateEF( StrandState& geometry, const ElasticStrand& strand )
{
    auto setIt = m_sets.find( &strand ) ;
    if ( setIt == m_sets.end() ) return ;

    const VerticesSet::IndexedSetType::SetType &strandVertices = setIt->second.getSet() ;

    Vec3x localF ;
    for( auto vIt = strandVertices.begin() ; vIt != strandVertices.end() ; ++vIt )
    {
       const Vec3x &v = geometry.getVertex( *vIt ) ;
       geometry.m_totalEnergy += computeLocalE( v, strand ) ;
       computeLocalF( v, strand, localF ) ;
       geometry.m_totalForce.segment< 3 >( 4 * *vIt ) += localF ;
    }
}

void PenaltyForce::accumulateE( StrandState& geometry, const ElasticStrand& strand )
{
    auto setIt = m_sets.find( &strand ) ;
    if ( setIt == m_sets.end() ) return ;

    const VerticesSet::IndexedSetType::SetType &strandVertices = setIt->second.getSet() ;

    for( auto vIt = strandVertices.begin() ; vIt != strandVertices.end() ; ++vIt )
    {
       const Vec3x &v = geometry.getVertex( *vIt ) ;
       geometry.m_totalEnergy += computeLocalE( v, strand ) ;
    }
}

void PenaltyForce::accumulateJ( StrandState& geometry, const ElasticStrand& strand )
{
    auto setIt = m_sets.find( &strand ) ;
    if ( setIt == m_sets.end() ) return ;

    const VerticesSet::IndexedSetType::SetType &strandVertices = setIt->second.getSet() ;

    Mat3x localJ ;
    for( auto vIt = strandVertices.begin() ; vIt != strandVertices.end() ; ++vIt )
    {
       const Vec3x &v = geometry.getVertex( *vIt ) ;
       computeLocalJ( v, strand, localJ ) ;
       geometry.m_totalJacobian->localStencilAdd< 3 >( 4 * *vIt, localJ ) ;
    }
}

void PenaltyForce::accumulateCurrentF( VecXx& globalF, const ElasticStrand& strand )
{
    auto setIt = m_sets.find( &strand ) ;
    if ( setIt == m_sets.end() ) return ;

    const VerticesSet::IndexedSetType::SetType &strandVertices = setIt->second.getSet() ;

    Vec3x localF ;
    for( auto vIt = strandVertices.begin() ; vIt != strandVertices.end() ; ++vIt )
    {
       const Vec3x &v = strand.getVertex( *vIt ) ;
       computeLocalF( v, strand, localF ) ;
       globalF.segment< 3 >( 4 * *vIt ) += localF ;
    }
}

void PenaltyForce::accumulateFutureE( Scalar& energy, const ElasticStrand& strand )
{
    int sidx = strand->getGlobalIndex();
    for( unsigned t = 0; t < m_bands->size(); ++t )
    {
        TwistEdge* edge = (*m_bands)[t];
        if( edge->parents.first->m_strand.getGlobalIndex() == sidx )
        {
            int v = edge->parents.first->m_vertexIndex;
            const Vec3x &vpos = strand.getVertex( v ) ;
            energy += computeLocalE( vpos, strand );
            const Vec3x &v1pos = strand.getVertex( v+1 ) ;
            energy += computeLocalE( v1pos, strand );
        }
        else if( edge->parents.second->m_strand.getGlobalIndex() == sidx )
        {

        }
    }
}

Scalar PenaltyForce::computeLocalE( const TwistEdge* edge, const ElasticStrand& strand )const
{
    const Scalar length = ( m_origin - v ).norm() ;
    return .5 * m_stiffness * square( length - m_thickness ) ;
}

void PenaltyForce::computeLocalF( const Vec3x &v, const ElasticStrand& strand, Vec3x& localF ) const
{
    const Scalar length = ( m_origin - v ).norm() ;
    if( !isSmall( length ) )
    {
        localF = m_stiffness * ( m_restLength / length - 1. ) * ( v - m_origin ) ;
    }
}

void PenaltyForce::computeLocalJ( const Vec3x &v, const ElasticStrand& strand, Mat3x &localJ ) const
{
   const Scalar length = ( m_origin - v ).norm() ;

   localJ.setIdentity() ;
   if( !isSmall( length ) )
   {
       const Vec3x e = ( m_origin - v ) / length ;
       localJ += m_restLength * ( e * e.transpose() - Mat3x::Identity() ) / length ;
   }

   localJ *= - m_stiffness ;
}
