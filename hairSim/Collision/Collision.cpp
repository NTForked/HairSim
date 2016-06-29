#include "Collision.h"
#include "VertexFaceCollision.h"
#include "EdgeFaceCollision.h"
#include "EdgeEdgeCollision.h"
#include "ElementProxy.h"
#include "CollisionUtils/CollisionUtils.h"
#include "../Strand/ElasticStrand.h"

static const double EXTRA_RADIUS = 1.e-3;

Collision::Collision( EdgeProxy* edgeProxy ):
    m_firstEdgeProxy( edgeProxy ) 
{}

Collision::~Collision()
{}

void Collision::print( std::ostream& os ) const
{
    os << "Print function not implemented yet";
}

std::ostream& operator<<( std::ostream& os, const Collision& col )
{
    col.print( os );
    return os;
}

bool compareTimes( const Collision* ct1, const Collision* ct2 )
{
    return ct1->m_time < ct2->m_time;
}

Vec3 FaceCollision::offset() const
{
    return ( m_offset.dot( m_normal ) + EXTRA_RADIUS ) * m_normal;
}

void CollidingPair::generateTransformationMatrix()
{
    Mat3x &E = m_transformationMatrix;
    E.col( 0 ) = m_normal;

    Vec3 t1( -m_normal[1], m_normal[0], 0 );
    const Scalar nt1 = t1.squaredNorm();

    if( isSmall( nt1 ) )
    {
        t1 = Vec3 ( 0, -m_normal[2], m_normal[1] );
        t1.normalize();
    } else {
        t1 /= std::sqrt( nt1 );
    }

    E.col(1) = t1;
    E.col(2) = m_normal.cross( t1 ).normalized();
}

void CollidingPair::swapIfNecessary()
{
    if ( objects.first.globalIndex == -1 || (
            objects.second.globalIndex != -1
            && objects.second.globalIndex < objects.first.globalIndex ) )
    {
        m_normal = -m_normal;
        std::swap( objects.first, objects.second );
    }
}
