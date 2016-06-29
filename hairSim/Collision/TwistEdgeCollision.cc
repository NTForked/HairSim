
#include "TwistEdgeCollision.hh"
// #include "../Core/ElasticStrand.hh"
// #include "ElementProxy.hh"
#include "CollisionUtils.hh"
#include "../Utils/Distances.hh"
#include "../Utils/TextLog.hh"
// #include "../Dynamic/StrandDynamics.hh"

// #include "../Render/OpenGLDecl.hh"
#include "../Dynamic/Config.hh"

#include "CTCD.h"
#include "CollisionDetector.hh"

static const double SQ_TOLERANCE = 1e-12;

using namespace std;

bool pruneCollision( teh )
{
    bool ignoreCollision = false;

    //temp, disabling collisions between bands
    if( m_firstProxy->isTwistedBand || m_secondProxy->isTwistedBand ){
        return false;
    }

    for( unsigned c = 0; c < m_firstProxy->children.size(); ++c )
    { // ignore if a band exists between these edges already
        if( m_firstProxy->children[c].first == m_secondProxy->uniqueID ) return false;
    }

    /// First check against reasons why we might outright ignore this contact
    if ( m_firstStrand == m_secondStrand && abs( m_firstVertex - m_secondVertex ) <= 1 ){
        return false; // dissallow collisions between vertices that form an edge. (directly neighboring Edges)
    }

    if( teh->m_frozenScene && m_firstProxy->frozenID != teh->m_frozenCheck && m_secondProxy->frozenID != teh->m_frozenCheck ){
        return false; // When rechecking, only test against those edges that have moved
    }

    if( rejectSelfCollisions && m_firstStrand->getGlobalIndex() == m_secondStrand->getGlobalIndex() ){
        return false; // config flag set to ignore self collisions...
    } 

    if( m_firstProxy->flagged && m_secondProxy->flagged ){
        return false;
    }

    if( m_secondProxy->isTwistedBand ){
        // ignore band collision with housing parents
        if( m_firstProxy == m_secondProxy->parents.first || m_firstProxy == m_secondProxy->parents.second ) return false;

        if( m_firstProxy == m_secondProxy->parents.first->prev || m_firstProxy == m_secondProxy->parents.first->next ) return false;
        if( m_firstProxy == m_secondProxy->parents.second->prev || m_firstProxy == m_secondProxy->parents.second->next ) return false;
    }

    // because not sorted have to do this twice
    if( m_firstProxy->isTwistedBand ){
        if( m_secondProxy == m_firstProxy->parents.first || m_secondProxy == m_firstProxy->parents.second ) return false;

        // ignore parent neighbors:
        if( m_secondProxy == m_firstProxy->parents.first->prev || m_secondProxy == m_firstProxy->parents.first->next ) return false;
        if( m_secondProxy == m_firstProxy->parents.second->prev || m_secondProxy == m_firstProxy->parents.second->next ) return false;

    }

    if( m_firstProxy->isTwistedBand && m_secondProxy->isTwistedBand ){
        return false;
    }


    return ignoreCollision;
}

bool TwistEdgeCollision::analyse( TwistEdgeHandler* teh )
{ // aka detectPotentiallyFrozenSceneEdgeEdgeNoThickness()

    if( pruneCollision( teh ) ){
        return false;
    }

    // post timestep positions
    Vec3 xA_final, xB_final;
    teh->getEdgeVerts( m_firstProxy, false, xA_final, xB_final );

    Vec3 xE_final, xF_final;
    teh->getEdgeVerts( m_secondProxy, false, xE_final, xF_final );

    // displacements (which give start positions)
    Vec3 dxA, dxB, dxE, dxF;
    if( teh->m_frozenScene ){
        if( m_firstProxy->frozenID == teh->m_frozenCheck ){
            Vec3 xA, xB; teh->getEdgeVerts( m_firstProxy, true, xA, xB );
            dxA = xA_final - xA;
            dxB = xB_final - xB;
        }
        else{
            dxA = Vec3::Zero();
            dxB = Vec3::Zero();
        }
        if( m_secondProxy->frozenID == teh->m_frozenCheck ){
            Vec3 xE, xF; teh->getEdgeVerts( m_secondProxy, true, xE, xF );
            dxE = xE_final - xE;
            dxF = xF_final - xF;
        }
        else{
            dxE = Vec3::Zero();
            dxF = Vec3::Zero();
        }
    }
    else{
        Vec3 xA, xB; teh->getEdgeVerts( m_firstProxy, true, xA, xB );
        dxA = xA_final - xA;
        dxB = xB_final - xB;
        Vec3 xE, xF; teh->getEdgeVerts( m_secondProxy, true, xE, xF );
        dxE = xE_final - xE;
        dxF = xF_final - xF;       
    }

    Vec3 cp, cq;
    double dist_squared = ClosestPtSegmentSegment( xA_final, xB_final, xE_final, xF_final, m_s, m_t, cp, cq );
    
    // Instantaneous Time Collision Detection
    Scalar radiusSquared_ofDesiredThickness = (m_firstProxy->m_radius + m_secondProxy->m_radius) * (m_firstProxy->m_radius + m_secondProxy->m_radius);
    if( dist_squared < radiusSquared_ofDesiredThickness )
    {
        m_time = 1.0;

        m_normal = cq - cp;// normalFunc( xA_final - dxA, xB_final - dxB, xE_final - dxE, xF_final - dxF, xA_final, xB_final, xE_final, xF_final, m_s, m_t );
        m_normal.normalize();

///////
// H // asdfasdfsfsdfsa Look at EdgeFaceIntersection code for Edge Edge proximity and see if that is of interest here....
        // also analyseRoughRodRodCollision
////////

        return true;
    }

    return false;
}
