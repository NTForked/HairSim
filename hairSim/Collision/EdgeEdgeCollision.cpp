#include "EdgeEdgeCollision.h"
#include "ElementProxy.h"
#include "CollisionUtils/CollisionUtils.h"
#include "../Math/Distances.hh"
#include "../Strand/ElasticStrand.h"
#include "../Strand/StrandDynamics.h"
#include "CollisionUtils/CTCD.h"

using namespace std;

bool EdgeEdgeCollision::validCollision( TwistEdgeHandler* teh )
{

    if( m_firstEdgeProxy->getStrandPointer()->getGlobalIndex() == m_secondEdgeProxy->getStrandPointer()->getGlobalIndex() 
            && abs( m_firstEdgeProxy->getVertexIndex() - m_secondEdgeProxy->getVertexIndex() ) <= 1 ){
        return false; // dissallow collisions between vertices that form an edge.
    }

    if( m_firstEdgeProxy->getStrandPointer()->collisionParameters().m_rejectSelfCollisions
            && m_firstEdgeProxy->getStrandPointer()->getGlobalIndex() == m_secondEdgeProxy->getStrandPointer()->getGlobalIndex() ){
        return false; // flag to ignore self collisions...
    } 


    // Bands:

    TwistEdge* firstProxy = dynamic_cast< TwistEdge* >( m_firstEdgeProxy );
    TwistEdge* secondProxy = dynamic_cast< TwistEdge* >( m_secondEdgeProxy );
    if( firstProxy == NULL || secondProxy == NULL ) return true;

    //temp, disabling collisions between bands
    if( firstProxy->isTwistedBand || secondProxy->isTwistedBand ){
        return false; // !!! Removing this requires switching get Start and End positions below to use teh->getEdgeVerts() func
    }

    for( unsigned c = 0; c < firstProxy->children.size(); ++c )
    { // ignore if a band exists between these edges already
        if( firstProxy->children[c].first == secondProxy->uniqueID ) return false;
    }

    if( teh->m_frozenScene && firstProxy->frozenID != teh->m_frozenCheck && secondProxy->frozenID != teh->m_frozenCheck ){
        return false; // When rechecking, only test against those edges that have moved
    }

    if( firstProxy->flagged && secondProxy->flagged ){
        return false;
    }

    if( secondProxy->isTwistedBand ){
        // ignore band collision with housing parents
        if( firstProxy == secondProxy->parents.first || firstProxy == secondProxy->parents.second ) return false;

        if( firstProxy == secondProxy->parents.first->prev || firstProxy == secondProxy->parents.first->next ) return false;
        if( firstProxy == secondProxy->parents.second->prev || firstProxy == secondProxy->parents.second->next ) return false;
    }

    // because not sorted have to do this twice
    if( firstProxy->isTwistedBand ){
        if( secondProxy == firstProxy->parents.first || secondProxy == firstProxy->parents.second ) return false;

        // ignore parent neighbors:
        if( secondProxy == firstProxy->parents.first->prev || secondProxy == firstProxy->parents.first->next ) return false;
        if( secondProxy == firstProxy->parents.second->prev || secondProxy == firstProxy->parents.second->next ) return false;

    }

    return true;
}

bool compare( const EdgeEdgeCollision* ef1, const EdgeEdgeCollision* ef2 )
{
    if( ef1->m_firstEdgeProxy == ef2->m_firstEdgeProxy ){
        return ef1->m_secondEdgeProxy < ef2->m_secondEdgeProxy;
    }
    else{
        return ef1->m_firstEdgeProxy < ef2->m_firstEdgeProxy;
    }
}

bool EdgeEdgeCollision::analyse( TwistEdgeHandler* teh )
{
    if( !validCollision( teh ) ){
        return false;
    }

    if( m_firstEdgeProxy->getStrandPointer()->collisionParameters().m_useInstantenousForCT && analyseInstantaneous() ){
        return true;
    }
    else if( !m_firstEdgeProxy->getStrandPointer()->collisionParameters().m_useInstantenousForCT && analyseContinousTime() ){
        return true;
    }

    return false;
}

/*
    Edges AB and CD, at start (s) and end (e) of timestep
*/
bool EdgeEdgeCollision::analyseContinousTime()
{
    // start of timestep positions
    const Vec3 pAs = m_firstEdgeProxy->getStrandPointer()->getVertex( m_firstEdgeProxy->getVertexIndex() );
    const Vec3 pBs = m_firstEdgeProxy->getStrandPointer()->getVertex( m_firstEdgeProxy->getVertexIndex() + 1 );
    const Vec3 pCs = m_secondEdgeProxy->getStrandPointer()->getVertex( m_secondEdgeProxy->getVertexIndex() );
    const Vec3 pDs = m_secondEdgeProxy->getStrandPointer()->getVertex( m_secondEdgeProxy->getVertexIndex() + 1 );

    // end of timestep positions
    const Vec3 pAe = m_firstEdgeProxy->getStrandPointer()->getFutureVertex( m_firstEdgeProxy->getVertexIndex() );
    const Vec3 pBe = m_firstEdgeProxy->getStrandPointer()->getFutureVertex( m_firstEdgeProxy->getVertexIndex() + 1 );
    const Vec3 pCe = m_secondEdgeProxy->getStrandPointer()->getFutureVertex( m_secondEdgeProxy->getVertexIndex() );
    const Vec3 pDe = m_secondEdgeProxy->getStrandPointer()->getFutureVertex( m_secondEdgeProxy->getVertexIndex() + 1 );

    Scalar desiredThickness = m_firstEdgeProxy->getStrandPointer()->collisionParameters().collisionRadius( m_firstEdgeProxy->getVertexIndex() ) 
                            + m_secondEdgeProxy->getStrandPointer()->collisionParameters().collisionRadius( m_secondEdgeProxy->getVertexIndex() );

    double colTime = -1.0;
    if( CTCD::edgeEdgeCTCD( pAs, pBs, pCs, pDs, pAe, pBe, pCe, pDe, desiredThickness, colTime ) ){

        // Time of collision
        const Vec3 pAc = (1.0 - colTime) * pAs + colTime * pAe;
        const Vec3 pBc = (1.0 - colTime) * pBs + colTime * pBe;
        const Vec3 pCc = (1.0 - colTime) * pCs + colTime * pCe;
        const Vec3 pDc = (1.0 - colTime) * pDs + colTime * pDe;

        Vec3 cAB, cCD;
        double dist_squared = ClosestPtSegmentSegment( pAc, pBc, pCc, pDc, m_s, m_t, cAB, cCD );

        m_time = colTime;
        m_normal = cCD - cAB;
        m_normal.normalize();

        return true;
    }

    return false;
}

bool EdgeEdgeCollision::analyseInstantaneous()
{
    // end of timestep positions
    const Vec3 pAe = m_firstEdgeProxy->getStrandPointer()->getFutureVertex( m_firstEdgeProxy->getVertexIndex() );
    const Vec3 pBe = m_firstEdgeProxy->getStrandPointer()->getFutureVertex( m_firstEdgeProxy->getVertexIndex() + 1 );
    const Vec3 pCe = m_secondEdgeProxy->getStrandPointer()->getFutureVertex( m_secondEdgeProxy->getVertexIndex() );
    const Vec3 pDe = m_secondEdgeProxy->getStrandPointer()->getFutureVertex( m_secondEdgeProxy->getVertexIndex() + 1 );

    Scalar desiredThickness = m_firstEdgeProxy->getStrandPointer()->collisionParameters().collisionRadius( m_firstEdgeProxy->getVertexIndex() ) 
                            + m_secondEdgeProxy->getStrandPointer()->collisionParameters().collisionRadius( m_secondEdgeProxy->getVertexIndex() );

    Vec3 cAB, cCD;
    double dist_squared = ClosestPtSegmentSegment( pAe, pBe, pCe, pDe, m_s, m_t, cAB, cCD );

    if( dist_squared < desiredThickness * desiredThickness ){
        m_time = 1.0;
        m_normal = cCD - cAB;
        m_normal.normalize();

        return true;
    }
    return false;
}

