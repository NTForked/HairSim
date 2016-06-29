#include "VertexFaceCollision.hh"
#include "../Core/ElasticStrand.hh"
#include "ElementProxy.hh"
#include "CollisionUtils.hh"
#include "../Utils/Distances.hh"
#include "../Utils/TextLog.hh"
#include "../Dynamic/StrandDynamics.hh"

static const double SQ_TOLERANCE = 1e-12;

void VertexFaceCollision::print( std::ostream& os ) const
{
    os << "VertexFaceCollision: strand vertex " << m_firstStrand->getGlobalIndex() << ' '
            << m_firstVertex << " vs. " << m_faceProxy << " by " << -m_normalRelativeDisplacement
            << " at time = " << m_time << '\n';
    os << "Normal displacement = " << m_normal << '\n';
    os << "Vertex moved from: " << m_firstStrand->getFutureState().getVertex( m_firstVertex )
            << " to " << m_firstStrand->getVertex( m_firstVertex ) << '\n';
    os << "Face moved: " << *m_faceProxy << '\n';
    os << "Normal relative displacement: " << m_normalRelativeDisplacement << '\n';
}

bool VertexFaceCollision::analyse()
{
    static __thread double times[4];

    const Vec3 p_off = m_firstStrand->getVertex( m_firstVertex );
    // NB after taking off the offset p = 0
    Vec3 pf0 = m_faceProxy->getVertex( 0 ) - p_off;
    Vec3 pf1 = m_faceProxy->getVertex( 1 ) - p_off;
    Vec3 pf2 = m_faceProxy->getVertex( 2 ) - p_off;

    const Vec3 dp = m_firstStrand->dynamics().getDisplacement( m_firstVertex );
    Vec3 df0 = m_faceProxy->getDisplacement( 0 );
    Vec3 df1 = m_faceProxy->getDisplacement( 1 );
    Vec3 df2 = m_faceProxy->getDisplacement( 2 );

    int fnsign = m_faceProxy->knowsNormalSign( true, *m_firstStrand, m_firstVertex ) ;
    m_normal = m_faceProxy->getNormal() ;


    const Scalar extraRadius = m_firstStrand->collisionParameters().externalCollisionsRadius( m_firstVertex, M_PI_2 ) ;

    if ( extraRadius < 0. )
        return false;

    const Vec3& prox = extraRadius * fnsign * m_normal;
    pf0 += prox;
    pf1 += prox;
    pf2 += prox;
    df0 += prox;
    df1 += prox;
    df2 += prox;

    unsigned num_times ;
    getCoplanarityTimes( -dp, pf0 - df0, pf1 - df1, pf2 - df2, Vec3(), pf0, pf1, pf2, times,
            NULL, num_times );

    for ( size_t j = 0; j < num_times; ++j )
    {
        // TODO: Use barycentric coordinates or point-triangle closest point < epsilon here? closest point < epsilon really just extends the triangle a little bit.
        // Determine if the collision actually happens
        const Scalar dtime = times[j] - 1.0;
        const Vec3 pcol = dtime * ( dp );
        const Vec3 f0col = pf0 + dtime * df0;
        const Vec3 f1col = pf1 + dtime * df1;
        const Vec3 f2col = pf2 + dtime * df2;

        const Vec3 cp = ClosestPtPointTriangle( pcol, f0col, f1col, f2col );

        // If, when they are coplanar, the objects are sufficiently close, register a collision
        if ( ( pcol - cp ).squaredNorm() < SQ_TOLERANCE )
        {
            m_time = times[j];

//            Scalar m_u, m_v, m_w;
            computeBarycentricCoordinates( f0col, f1col, f2col, pcol, m_u, m_v, m_w );
            // computeBarycentricCoordinates coords could be outside of [0,1] right now because we've extended the triangles a little bit
            assert( isSmall(m_u + m_v + m_w - 1.0) );

            m_meshDisplacement = m_u * df0 + m_v * df1 + m_w * df2;
            const Vec3 relativeDisplacement = ( 1 - m_time ) * ( dp - m_meshDisplacement );

            //            m_offset = relativeDisplacement ;
            m_offset = ( m_u * pf0 + m_v * pf1 + m_w * pf2 ) - m_meshDisplacement // orig point on mesh
                    + dp;  // orig point on rod


            const Scalar nDepl = relativeDisplacement.dot( m_normal ) ;
            if( !fnsign )
            {
                // Normal sign was unknown, now we know that it should be opposed to relativeDisplacement
                fnsign = ( nDepl > 0. ? -1 : 1 ) ;
                m_faceProxy->setNormalSign( fnsign, m_time, *m_firstStrand, m_firstVertex );
//                m_meshDisplacement += extraRadius * fnsign * m_normal ;
            }
            else {
                if ( fnsign * nDepl > 0. )
                {
                    return false;
                }

                const Scalar belowRadius = m_offset.dot( fnsign * m_normal ) + extraRadius ;
                if( belowRadius > 0. && extraRadius > 0. )
                {
                    // We're already closer than the collision radius, let proximity handle that
                    m_meshDisplacement -= prox ;
                    m_offset.setZero() ;
                }

            }
            m_normal = fnsign * m_normal ;


            postAnalyse( relativeDisplacement );

            return true;
        }
    }
    return false;

}

bool compare( const VertexFaceCollision* vf1, const VertexFaceCollision* vf2 )
{
    if ( vf1->m_firstStrand == vf2->m_firstStrand )
        if ( vf1->m_firstVertex == vf2->m_firstVertex )
            if ( vf1->m_faceProxy == vf2->m_faceProxy )
                return vf1->m_time < vf2->m_time;
            else
                return vf1->m_faceProxy < vf2->m_faceProxy;
        else
            return vf1->m_firstVertex < vf2->m_firstVertex;
    else
        return vf1->m_firstStrand < vf2->m_firstStrand;
}
