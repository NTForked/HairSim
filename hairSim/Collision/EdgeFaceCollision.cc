#include "EdgeFaceCollision.hh"
#include "../Core/ElasticStrand.hh"
#include "ElementProxy.hh"
#include "CollisionUtils.hh"
#include "../Utils/Distances.hh"
#include "../Utils/TextLog.hh"
#include "../Dynamic/StrandDynamicTraits.hh"


static const double SQ_TOLERANCE = 1e-12;

void EdgeFaceCollision::print( std::ostream& os ) const
{
    os << "EdgeFaceCollision: strand edge " << m_firstStrand->getGlobalIndex() << ' ' << m_firstVertex
            << " s = " << m_s;
    os << " vs. face edge " << m_mesh << ' ' << m_firstIdx << ' ' << m_secondIdx << " by "
            << -m_normalRelativeDisplacement << " at time = " << m_time;

    // Extra debugging info if needed
    os << '\n';
    os << "Normal displacement = " << m_normal << '\n';

    os << "Strand edge moved from: " << m_firstStrand->getFutureState().getVertex( m_firstVertex )
            << " --- " << m_firstStrand->getFutureState().getVertex( m_firstVertex + 1 ) << ' ';
    os << "to " << m_firstStrand->getVertex( m_firstVertex ) << " --- "
            << m_firstStrand->getVertex( m_firstVertex + 1 );
    os << '\n';

    os << "Face edge moved from: "
            << m_mesh->getVertex( m_firstIdx ) - m_mesh->getDisplacement( m_firstIdx ) << " --- "
            << m_mesh->getVertex( m_secondIdx ) - m_mesh->getDisplacement( m_secondIdx ) << ' ';
    os << "to " << m_mesh->getVertex( m_firstIdx ) << " --- " << m_mesh->getVertex( m_secondIdx );
    os << '\n';

    os << "Normal relative displacement: " << m_normalRelativeDisplacement;
    os << '\n';
}

bool EdgeFaceCollision::analyse()
{
    static __thread double times[4];
//    const short side_p = m_side < 2 ? m_side + 1 : 0;

    // Tolerance for the interior edge/edge collisions that point in a direction perpendicular to the face normal
    const Scalar perpEdgeTol = .1 ;

    const Vec3x pp0 = m_firstStrand->getVertex( m_firstVertex );
    const Vec3x pp1 = m_firstStrand->getVertex( m_firstVertex + 1 );
    Vec3x pq0 = m_mesh->getVertex( m_firstIdx );
    Vec3x pq1 = m_mesh->getVertex( m_secondIdx );

    Vec3x dp0 = m_firstStrand->dynamics().getDisplacement( m_firstVertex );
    Vec3x dp1 = m_firstStrand->dynamics().getDisplacement( m_firstVertex + 1 );
    Vec3x dq0 = m_mesh->getDisplacement( m_firstIdx );
    Vec3x dq1 = m_mesh->getDisplacement( m_secondIdx );

    int fnsign = 0 ;
    Vec3x faceNormal ;

    bool shouldAddProx = false ;

    if( m_onBoundary )
    {
        const short thirdApex = 3 - m_firstApex - m_secondApex ;
        const Vec3x pq2 = m_mesh->getVertex(  m_face->getFace().idx[ thirdApex ] );

        const Vec3x meshEdge = ( pq1 - pq0 ).normalized() ;
        faceNormal = ( pq0 - pq2 ) - ( pq0 - pq2 ).dot( meshEdge ) * meshEdge ;

        const Scalar nnorm = faceNormal.norm();
        if( isSmall( nnorm ) ) return false ;
        faceNormal /= nnorm ;

    } else {

        fnsign = m_face->knowsNormalSign( true, *m_firstStrand, m_firstVertex ) ;
        faceNormal = m_face->getNormal() ;

        // Try to guess approximate collision normal so we can try to enforce collision radius
        m_normal = ( pq1 - pq0 ).cross( ( pp1 - dp1 ) - ( pp0 - dp1 ) );
        const Scalar nnorm = m_normal.norm();

        if( fnsign )
        {
            shouldAddProx = true ;

            if ( isSmall( nnorm ) )
            {
                m_normal = faceNormal * fnsign  ;
            }
            else
            {
                m_normal /= nnorm ;
                if( fnsign * m_normal.dot( faceNormal ) < 0. )
                {
                    m_normal *= -1.;
                }
            }
        } else {

            shouldAddProx = ! isSmall( nnorm ) ;
            if ( shouldAddProx )
            {
                m_normal /= nnorm ;

                Scalar s = .5, t = .5 ;
                const Vec3x meshDisplacement = ( 1.0 - t ) * dq0 + t * dq1;
                const Vec3x relativeDisplacement =
                        ( ( ( 1.0 - s ) * dp0 + s * dp1 ) - meshDisplacement );
                if( m_normal.dot( relativeDisplacement ) > 0. )
                {
                    m_normal *= -1. ;
                }

            }
        }
        shouldAddProx = shouldAddProx && std::fabs( faceNormal.dot( m_normal ) ) > perpEdgeTol ;
    }


    const Scalar extraRadius = m_firstStrand->collisionParameters().externalCollisionsRadius( m_firstVertex, M_PI_2 ) ;

    if ( extraRadius < 0. )
        return false;

    Vec3x prox ;
    if( shouldAddProx )
    {
        prox = extraRadius * m_normal;
        pq0 += prox;
        pq1 += prox;
        dq0 += prox;
        dq1 += prox;
    }

    // TODO: tetrahedron test?

    unsigned num_times ;
    getCoplanarityTimes( pp0 - dp0, pp1 - dp1, pq0 - dq0, pq1 - dq1, pp0, pp1, pq0, pq1, times,
            NULL, num_times );

    // Loop over the coplanarity times until we find a bona fide collision
    for ( unsigned j = 0 ; j < num_times ; ++j )
    {
        const Scalar dtime = times[j] - 1.0;

        // Determine if the collision actually happens
        const Vec3x p0col = pp0 + dtime * ( dp0 );
        const Vec3x p1col = pp1 + dtime * ( dp1 );
        const Vec3x q0col = pq0 + dtime * dq0;
        const Vec3x q1col = pq1 + dtime * dq1;

        Vec3x cp, cq;
        Scalar t;
        const double sqrdist = ClosestPtSegmentSegment( p0col, p1col, q0col, q1col, m_s, t, cp,
                cq );

        // Compute the barycentric coordinates of the hit, even though it's on the side of the triangle
        m_u = m_v = m_w = 0;
        switch ( m_firstApex )
        {
        case 0:
            m_u = 1 - t;
            break;
        case 1:
            m_v = 1 - t;
            break;
        case 2:
            m_w = 1 - t;
            break;
        }
        switch ( m_secondApex )
        {
        case 0:
            m_u = t;
            break;
        case 1:
            m_v = t;
            break;
        case 2:
            m_w = t;
            break;
        }

        // If, when they are coplanar, the objects are sufficiently close, register a collision
        if ( sqrdist < SQ_TOLERANCE )
        {
            m_time = times[j];


            // Compute a collision normal at the time of the collision. For a first attempt, take the cross product of the edges.
//            m_normal = ( q1col - q0col ).cross( p1col - p0col );
            m_normal = ( pq1 - pq0 ).cross( pp1 - pp0 );
            double nnorm = m_normal.norm();


            // If the edges happen to be parallel
            if ( nnorm * nnorm <= SQ_TOLERANCE )
            {
                TraceStream( g_log, "" ) << "Parallel edges in EdgeFaceCollision";
                // Use the pre-timestep positions of the collision points to generate a collision normal
                m_normal = ( ( 1.0 - t ) * pq0 + t * pq1 ).cross( ( 1.0 - m_s ) * pp0 + m_s * pp1 );
                nnorm = m_normal.norm();

                assert( nnorm*nnorm > SQ_TOLERANCE );
            }
            m_normal /= nnorm;

            m_meshDisplacement = ( 1.0 - t ) * dq0 + t * dq1;

            const Vec3x relativeDisplacement = ( 1.0 - m_time )
                    * ( ( ( 1.0 - m_s ) * dp0 + m_s * dp1 ) - m_meshDisplacement );
            postAnalyse( relativeDisplacement );


            m_offset = ( ( 1.0 - t ) * pq0 + t * pq1 ) - m_meshDisplacement //p mesh orig
                    -  ( ( 1.0 - m_s ) * ( pp0 - dp0 ) + m_s * ( pp1 -dp1 ) ) ; //p rod orig


            Scalar perp = m_normal.dot( faceNormal ) ;

            if ( m_onBoundary )
            {
//                std::cout << m_firstVertex << " / " << perp << " / " << faceNormal << " " << m_normal << std::endl ;
                if( perp < 0 )
                {
                    if( perp < -.9 ) continue ; // Collision normal opposed to boundary normal, discard

                    m_normal -= m_normal.dot( faceNormal ) * faceNormal ;
                    m_normal = m_normal.normalized() ;
                }
//                m_meshDisplacement += extraRadius * m_normal ;

            } else if( fnsign ) {
                if( perp * fnsign < 0. )
                {
                    m_normal -= 2 * perp * faceNormal ;
                    perp = -perp ;
                }

                // Collision on an inside edge almost perpendicular to face normal, discard
                if( perp * fnsign < perpEdgeTol )
                    continue ;
            }
            else
            {
                // Collision against an inner edge, we can deduce face normal
                const Scalar perpDepl = relativeDisplacement.dot( faceNormal) ;

                if( perp * perpDepl > 0 ) continue ; // Relative displacememnt and collison normal disagree, ignore sign

                fnsign = ( perpDepl < 0. ? 1 : -1 ) ;

                m_face->setNormalSign( fnsign, m_time, *m_firstStrand, m_firstVertex );

                //                m_meshDisplacement += extraRadius * m_normal ;

            }

            if( shouldAddProx )
            {
                const Scalar extraRadius = prox.dot( m_normal ) ;
                const Scalar belowRadius = m_offset.dot( m_normal ) + extraRadius ;
                if( belowRadius > 0. && extraRadius > 0. )
                {
                    // We're already closer than the collision radius, let proximity handle that
                    m_meshDisplacement -= prox ;
                    m_offset.setZero() ;
                }
            }

            return true;
        }
    }
    return false;

}

bool compare( const EdgeFaceCollision* ef1, const EdgeFaceCollision* ef2 )
{
    if ( ef1->m_firstStrand == ef2->m_firstStrand )
        if ( ef1->m_firstVertex == ef2->m_firstVertex )
            if ( ef1->m_firstIdx == ef2->m_firstIdx )
                if ( ef1->m_secondIdx == ef2->m_secondIdx )
                    return ef1->m_time < ef2->m_time;
                else
                    return ef1->m_secondIdx < ef2->m_secondIdx;
            else
                return ef1->m_firstIdx < ef2->m_firstIdx;
        else
            return ef1->m_firstVertex < ef2->m_firstVertex;
    else
        return ef1->m_firstStrand < ef2->m_firstStrand;
}
