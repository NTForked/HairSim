#include "EdgeEdgeCollision.hh"
#include "../Core/ElasticStrand.hh"
#include "ElementProxy.hh"
#include "CollisionUtils.hh"
#include "../Utils/Distances.hh"
#include "../Utils/TextLog.hh"
#include "../Dynamic/StrandDynamicTraits.hh"
#include "../Render/OpenGLDecl.hh"
#include "../Dynamic/Config.hh"
#include "CTCD.h"
#include "CollisionDetector.hh"

using namespace std;

bool EdgeEdgeCollision::analyse()
{
    if( m_firstStrand == m_secondStrand && abs( m_firstVertex - m_secondVertex ) <= 1 ){
        return false; // dissallow collisions between vertices that form an edge.
    }

    // post timestep positions
    const Vec3x pp0 = m_firstStrand->getVertex( m_firstVertex );
    const Vec3x pp1 = m_firstStrand->getVertex( m_firstVertex + 1 );
    const Vec3x pq0 = m_secondStrand->getVertex( m_secondVertex );
    const Vec3x pq1 = m_secondStrand->getVertex( m_secondVertex + 1 );

    // displacements
    const Vec3x dp0 = m_firstStrand->dynamics().getDisplacement( m_firstVertex );
    const Vec3x dp1 = m_firstStrand->dynamics().getDisplacement( m_firstVertex + 1 );
    const Vec3x dq0 = m_secondStrand->dynamics().getDisplacement( m_secondVertex );
    const Vec3x dq1 = m_secondStrand->dynamics().getDisplacement( m_secondVertex + 1 );


    // Continous Time collision Detection
    double times[4];
    unsigned num_times ;
    getCoplanarityTimes( pp0 - dp0, pp1 - dp1, pq0 - dq0, pq1 - dq1, pp0, pp1, pq0, pq1, times, NULL, num_times );

    // Loop over the coplanarity times until we find a bona fide collision
    for ( unsigned j = 0; j < num_times ; ++j )
    {
        const Scalar dtime = times[j] - 1.0;

        // Determine if the collision actually happens
        const Vec3x p0col = pp0 + dtime * dp0;
        const Vec3x p1col = pp1 + dtime * dp1;
        const Vec3x q0col = pq0 + dtime * dq0;
        const Vec3x q1col = pq1 + dtime * dq1;

        Vec3x cp, cq;
        double dist_squared = ClosestPtSegmentSegment( p0col, p1col, q0col, q1col, m_s, m_t, cp, cq );
        
        Scalar radiusSquared_ofDesiredThickness = SQ_TOLERANCE;
        if( dist_squared < radiusSquared_ofDesiredThickness )
        {

            if( rejectSelfCollisions && m_firstStrand->getGlobalIndex() == m_secondStrand->getGlobalIndex() ){
                return false; // Temp ignoring other collisions...
            }

            m_time = times[j];


            std::cout << "Need to fix normal here: " << m_time <<endl;


            m_normal = cq - cp;// normalFunc( pp0 - dp0, pp1 - dp1, pq0 - dq0, pq1 - dq1, p0col, p1col, q0col, q1col, m_s, m_t );
            double nnorm = m_normal.norm();

            m_normal.normalize();

            // vector pointing from one collision point to other
            const Vec3x relativeDisplacement = ( 1.0 - m_time ) * ( ( ( 1.0 - m_s ) * dp0 + m_s * dp1 ) - ( ( 1.0 - m_t ) * dq0 + m_t * dq1 ) );
            postAnalyse( relativeDisplacement );

            m_penetrationVector = Vec3x::Zero();

            if( thickness ){
                m_normal =  -m_normal; // usually necesary when thickness > 1e-3....
            }

            if( !velCondition && !penetrationResponse ){
                return true;
            }

            Vec3x x10 = ((1 - m_s) * (pp0 - dp0) + m_s * (pp1 - dp1) );
            Vec3x x20 = ((1 - m_t) * (pq0 - dq0) + m_t * (pq1 - dq1) );
            const Scalar normal_dist0   = (x20 - x10).dot(m_normal);
            const Scalar min_normal_vel = -(sqrt(radiusSquared_ofDesiredThickness ) - normal_dist0) / 1e-3;

            if( penetrationResponse ){
                m_penetrationVector = min_normal_vel * 1e-3 * 1e-3 * m_normal; // NEED timestep m_dt here
            }

            if( !velCondition ){
                return true;        
            }
            const Scalar residual_velocity = 1e-6;

            const Vec3x vp0 = m_firstStepper->newVelocities().segment<3>(    m_firstVertex       * 4 );
            const Vec3x vp1 = m_firstStepper->newVelocities().segment<3>(   (m_firstVertex + 1 ) * 4 );
            const Vec3x vq0 = m_secondStepper->newVelocities().segment<3>(  m_secondVertex       * 4 );
            const Vec3x vq1 = m_secondStepper->newVelocities().segment<3>( (m_secondVertex + 1 ) * 4 );

            Vec3x v1  = ((1 - m_s) * vp0 + m_s * vp1);
            Vec3x v2  = ((1 - m_t) * vq0 + m_t * vq1);

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
