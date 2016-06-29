#include "EdgeEdgeCollision.h"
#include "ElementProxy.h"
#include "CollisionUtils/CollisionUtils.h"
#include "../Math/Distances.hh"
#include "../Strand/ElasticStrand.h"
#include "../Strand/StrandDynamics.h"

using namespace std;

bool EdgeEdgeCollision::validCollision( TwistEdgeHandler* teh )
{

    if( m_firstStrand == m_secondStrand && abs( m_firstVertex - m_secondVertex ) <= 1 ){
        return false; // dissallow collisions between vertices that form an edge.
    }

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

    if( m_firstStrand.params.rejectSelfCollisions && m_firstStrand->getGlobalIndex() == m_secondStrand->getGlobalIndex() ){
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

    return true;
}



bool EdgeEdgeCollision::analyse( TwistEdgeHandler* teh  )
{
    if( !validCollision( teh ) ){
        return false;
    }

    // post timestep positions
    const Vec3 pp0 = m_firstStrand->getVertex( m_firstVertex );
    const Vec3 pp1 = m_firstStrand->getVertex( m_firstVertex + 1 );
    const Vec3 pq0 = m_secondStrand->getVertex( m_secondVertex );
    const Vec3 pq1 = m_secondStrand->getVertex( m_secondVertex + 1 );

    // displacements
    const Vec3 dp0 = m_firstStrand->dynamics().getDisplacement( m_firstVertex );
    const Vec3 dp1 = m_firstStrand->dynamics().getDisplacement( m_firstVertex + 1 );
    const Vec3 dq0 = m_secondStrand->dynamics().getDisplacement( m_secondVertex );
    const Vec3 dq1 = m_secondStrand->dynamics().getDisplacement( m_secondVertex + 1 );

    // Continous Time collision Detection
    double times[4];
    unsigned num_times ;
    getCoplanarityTimes( pp0 - dp0, pp1 - dp1, pq0 - dq0, pq1 - dq1, pp0, pp1, pq0, pq1, times, NULL, num_times );

    // Loop over the coplanarity times until we find a bona fide collision
    for ( unsigned j = 0; j < num_times ; ++j )
    {
        const Scalar dtime = times[j] - 1.0;

        // Determine if the collision actually happens
        const Vec3 p0col = pp0 + dtime * dp0;
        const Vec3 p1col = pp1 + dtime * dp1;
        const Vec3 q0col = pq0 + dtime * dq0;
        const Vec3 q1col = pq1 + dtime * dq1;

        Vec3 cp, cq;
        double dist_squared = ClosestPtSegmentSegment( p0col, p1col, q0col, q1col, m_s, m_t, cp, cq );
        
        Scalar radiusSquared_ofDesiredThickness = strand.cparams.collisionRadius(m_firstVertex) + strand.cparams.collisionRadius( m_secondVertex );
        radiusSquared_ofDesiredThickness *= radiusSquared_ofDesiredThickness;
        if( dist_squared < radiusSquared_ofDesiredThickness )
        {

            if( rejectSelfCollisions && m_firstStrand->getGlobalIndex() == m_secondStrand->getGlobalIndex() ){
                return false; // Temp ignoring other collisions...
            }

            m_time = times[j];

            m_normal = cq - cp;
            m_normal.normalize();

            // vector pointing from one collision point to other
            const Vec3 relativeDisplacement = ( 1.0 - m_time ) * ( ( ( 1.0 - m_s ) * dp0 + m_s * dp1 ) - ( ( 1.0 - m_t ) * dq0 + m_t * dq1 ) );
            postAnalyse( relativeDisplacement );

            m_penetrationVector = Vec3::Zero();

            if( thickness ){
                m_normal =  -m_normal; // usually necesary when thickness > 1e-3....
            }

            if( !velCondition && !penetrationResponse ){
                return true;
            }

            Vec3 x10 = ((1 - m_s) * (pp0 - dp0) + m_s * (pp1 - dp1) );
            Vec3 x20 = ((1 - m_t) * (pq0 - dq0) + m_t * (pq1 - dq1) );
            const Scalar normal_dist0   = (x20 - x10).dot(m_normal);
            const Scalar min_normal_vel = -(sqrt(radiusSquared_ofDesiredThickness ) - normal_dist0) / 1e-3;

            if( penetrationResponse ){
                m_penetrationVector = min_normal_vel * 1e-3 * 1e-3 * m_normal; // NEED timestep m_dt here
            }

            if( !velCondition ){
                return true;        
            }
            const Scalar residual_velocity = 1e-6;

            const Vec3 vp0 = m_firstStepper->newVelocities().segment<3>(    m_firstVertex       * 4 );
            const Vec3 vp1 = m_firstStepper->newVelocities().segment<3>(   (m_firstVertex + 1 ) * 4 );
            const Vec3 vq0 = m_secondStepper->newVelocities().segment<3>(  m_secondVertex       * 4 );
            const Vec3 vq1 = m_secondStepper->newVelocities().segment<3>( (m_secondVertex + 1 ) * 4 );

            Vec3 v1  = ((1 - m_s) * vp0 + m_s * vp1);
            Vec3 v2  = ((1 - m_t) * vq0 + m_t * vq1);

            Scalar normal_vel;
            if( thickness ){
                normal_vel = (v2  - v1 ).dot(m_normal);
            }
            else{
                normal_vel = (v1  - v2 ).dot(m_normal); // Usually necessary when SQ_TOLERANCE 
            }

#pragma omp critical
    {       
            if (  normal_vel >= ( min_normal_vel - residual_velocity ) ) {    
                std::cerr << "EdgeEdgeCollision::analyse DISCARDING EEC " <<
                 m_firstStrand->getGlobalIndex() <<" " << m_secondStrand->getGlobalIndex() <<
                 ", moving away with normal_vel: " << normal_vel  << " which is  >= min_normal_vel: " << min_normal_vel << std::endl;
            }
            else{
                std::cerr << "EEC " << (  normal_vel < ( min_normal_vel - residual_velocity ) ) << ", moving with normal_vel: " << normal_vel << std::endl;
            }     
    }

            // If we are moving away, then do not return as registered collision.
            return (  normal_vel < ( min_normal_vel - residual_velocity ) );

        }
    }

    return false;
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

void EdgeEdgeCollision::print( std::ostream& os ) const
{
    os << "EdgeEdgeCollision: strand edge " << m_firstStrand->getGlobalIndex() << ' ' << m_firstVertex
        << " vs. strand edge " << m_secondStrand->getGlobalIndex() << ' ' << m_secondVertex
        << " by " << -m_normalRelativeDisplacement << " at time = " << m_time << " normal: " << m_normal; 
}

void EdgeEdgeCollision::printShort( std::ostream& os ) const
{
    os << "[" << m_firstStrand->getGlobalIndex() << "." << m_firstVertex << "|" <<
        m_secondStrand->getGlobalIndex() << "." << m_secondVertex<< "] time: " << m_time 
        << " N: " << m_normal;
}
