
#include "TwistEdgeCollision.hh"
// #include "../Core/ElasticStrand.hh"
// #include "ElementProxy.hh"
#include "CollisionUtils.hh"
#include "../Utils/Distances.hh"
#include "../Utils/TextLog.hh"
// #include "../Dynamic/StrandDynamicTraits.hh"

// #include "../Render/OpenGLDecl.hh"
#include "../Dynamic/Config.hh"

#include "CTCD.h"
#include "CollisionDetector.hh"

static const double SQ_TOLERANCE = 1e-12;

using namespace std;
namespace strandsim
{

bool TwistEdgeCollision::analyse( TwistEdgeHandler* teh )
{ // aka detectPotentiallyFrozenSceneEdgeEdgeNoThickness()

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
// cout << "frozen check " << endl;
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
// cout << "m_secondProxy parents first: " << m_secondProxy->parents.first->uniqueID << " second: " << m_secondProxy->parents.second->uniqueID << endl;
        if( m_firstProxy == m_secondProxy->parents.first || m_firstProxy == m_secondProxy->parents.second ) return false;

        if( m_firstProxy == m_secondProxy->parents.first->prev || m_firstProxy == m_secondProxy->parents.first->next ) return false;
        if( m_firstProxy == m_secondProxy->parents.second->prev || m_firstProxy == m_secondProxy->parents.second->next ) return false;
    }

    // because not sorted have to do this twice
    if( m_firstProxy->isTwistedBand ){
// cout << "m_firstProxy parents first: " << m_firstProxy->parents.first->uniqueID << " second: " << m_firstProxy->parents.second->uniqueID << endl;
        if( m_secondProxy == m_firstProxy->parents.first || m_secondProxy == m_firstProxy->parents.second ) return false;

        // ignore parent neighbors:
        if( m_secondProxy == m_firstProxy->parents.first->prev || m_secondProxy == m_firstProxy->parents.first->next ) return false;
        if( m_secondProxy == m_firstProxy->parents.second->prev || m_secondProxy == m_firstProxy->parents.second->next ) return false;

    }

    if( m_firstProxy->isTwistedBand && m_secondProxy->isTwistedBand ){
        return false;
    }


    // post timestep positions
    Vec3x xA_final, xB_final;
    teh->getEdgeVerts( m_firstProxy, false, xA_final, xB_final );
    Vec3x xE_final, xF_final;
    teh->getEdgeVerts( m_secondProxy, false, xE_final, xF_final );

    // displacements (which give start positions)
    Vec3x dxA, dxB, dxE, dxF;
    if( teh->m_frozenScene ){
        if( m_firstProxy->frozenID == teh->m_frozenCheck ){
            Vec3x xA, xB; teh->getEdgeVerts( m_firstProxy, true, xA, xB );
            dxA = xA_final - xA;
            dxB = xB_final - xB;
        }
        else{
            dxA = Vec3x::Zero();
            dxB = Vec3x::Zero();
        }
        if( m_secondProxy->frozenID == teh->m_frozenCheck ){
            Vec3x xE, xF; teh->getEdgeVerts( m_secondProxy, true, xE, xF );
            dxE = xE_final - xE;
            dxF = xF_final - xF;
        }
        else{
            dxE = Vec3x::Zero();
            dxF = Vec3x::Zero();
        }
    }
    else{
        Vec3x xA, xB; teh->getEdgeVerts( m_firstProxy, true, xA, xB );
        dxA = xA_final - xA;
        dxB = xB_final - xB;
        Vec3x xE, xF; teh->getEdgeVerts( m_secondProxy, true, xE, xF );
        dxE = xE_final - xE;
        dxF = xF_final - xF;       
    }

    Vec3x cp, cq;
    double dist_squared = ClosestPtSegmentSegment( xA_final, xB_final, xE_final, xF_final, m_s, m_t, cp, cq );
    
    // Scalar radiusSquared_ofDesiredThickness = SQ_TOLERANCE;
    Scalar radiusSquared_ofDesiredThickness = (m_firstProxy->radius + m_secondProxy->radius) * (m_firstProxy->radius + m_secondProxy->radius);
    // std::cout << "radiusSquared_ofDesiredThickness: " << radiusSquared_ofDesiredThickness << std::endl;
    if( dist_squared < radiusSquared_ofDesiredThickness )
    {
        m_time = 1.0;
    // std::cout << "radiusSquared_ofDesiredThickness: " << radiusSquared_ofDesiredThickness << " dist_squared: " << dist_squared << std::endl;

        m_normal = normalFunc( xA_final - dxA, xB_final - dxB, xE_final - dxE, xF_final - dxF, xA_final, xB_final, xE_final, xF_final, m_s, m_t );
        double nnorm = m_normal.norm();

        // If the edges happen to be parallel
        if ( nnorm * nnorm <= SQ_TOLERANCE )
        {
#pragma omp critical
            {
                std::cout << "ERROR, parallel collision picked up from edgeEdge old_analyse()" << std::endl;
                std::cout << "EdgeEdgeCollision: strand edge " << m_firstStrand->getGlobalIndex() << ' ' << m_firstVertex
                << " vs. strand edge " << m_secondStrand->getGlobalIndex() << ' ' << m_secondVertex
                << " by " << -m_normalRelativeDisplacement << " at time = " << m_time << " normal: " << m_normal << std::endl; 
            }
            // The comment on pre-timestep is WRONG. its using post-timestep
            // also it is taking the cross product of two points
            std::exit ( EXIT_FAILURE );
        }




///////
// H // asdfasdfsfsdfsa Look at EdgeFaceIntersection code for Edge Edge proximity and see if that is of interest here....
////////


        m_normal.normalize();
        m_normal *= -1.0;

        const Vec3x relativeDisplacement = ( 1.0 - m_time ) * ( ( ( 1.0 - m_s ) * dxA + m_s * dxB ) - ( ( 1.0 - m_t ) * dxE + m_t * dxF ) );
        postAnalyse( relativeDisplacement );

        m_penetrationVector = Vec3x::Zero();
        if( !penetrationResponse ){
            return true;
        }

        Vec3x x10 = ((1 - m_s) * (xA_final - dxA) + m_s * (xB_final - dxB) );
        Vec3x x20 = ((1 - m_t) * (xE_final - dxE) + m_t * (xF_final - dxF) );
        const Scalar normal_dist0   = (x20 - x10).dot(m_normal);
        const Scalar min_normal_vel = (sqrt(radiusSquared_ofDesiredThickness ) - normal_dist0);
        // cout << "should the calc above be + or - the thickness and the normal_dist0? " << endl;
        if( penetrationResponse ){
            m_penetrationVector = min_normal_vel * m_normal; // NEED timestep m_dt^2 here
        }

        if( m_penetrationVector.dot(m_normal) < 0.0 ){
            std::cout << "should be positive: collision pen dot normal " << m_penetrationVector.dot(m_normal) << std::endl;
            std::exit ( EXIT_FAILURE );
        }

        return true;
    }

    return false;
}

} // namespace strandsim
