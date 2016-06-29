#include "Simulation.h"

void Simulation::deleteInvertedProxies()
{
    if( penaltyAfter )
    {
        for ( std::vector<ElasticStrand*>::size_type i = 0; i < m_strands.size(); i++ )
        {
            m_collisionDetector->m_proxyHistory->applyImpulses( m_strands[i], m_steppers[i], true );
            m_steppers[i]->update( true );
        }
    }

    // grab just the tunneled bands, none of the originals
    std::vector< TwistEdge* > tunneledBands = m_collisionDetector->m_proxyHistory->tunnelingBands;
    for( unsigned t = 0; t < tunneledBands.size(); ++t )
    {

        TwistEdge* edge = tunneledBands[t];

        if( penaltyOnce ){
            m_collisionDetector->m_proxyHistory->deleteBand( edge, m_collisionDetector->m_elementProxies );                
            continue;
        }

        if( edge->traversed ){
            edge->traversed = false;
        }
        else{
            m_collisionDetector->m_proxyHistory->updateTwistAngle( edge, edge->parents.first, edge->parents.second, edge->parents.first, edge->parents.second, false );
        }

        if( edge->intersectionTwists() == 0 )
        {
            //check thickness boundary...
            Vec3 edgeA, edgeB;
            m_collisionDetector->m_proxyHistory->getEdgeVerts( edge, false, edgeA, edgeB );

            if( (edgeB - edgeA).norm() >= (2.0 * edge->m_radius) ){
                m_collisionDetector->m_proxyHistory->deleteBand( edge, m_collisionDetector->m_elementProxies );                
            }
        }
    }  
}
