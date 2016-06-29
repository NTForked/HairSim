#include "ContinuousTimeCollision.hh"
#include "VertexFaceCollision.hh"
#include "EdgeFaceCollision.hh"
#include "EdgeEdgeCollision.hh"
#include "ElementProxy.hh"
#include "CollisionUtils.hh"
#include "../Core/ElasticStrand.hh"

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

Vec3x FaceCollision::offset() const
{
    return ( m_offset.dot( m_normal ) + EXTRA_RADIUS ) * m_normal;
}



#include "ProximityCollision.hh"

void ProximityCollision::generateTransformationMatrix()
{
    Mat3x &E = m_transformationMatrix;
    E.col( 0 ) = normal ;

    Vec3x t1( -normal[1], normal[0], 0 ) ;
    const Scalar nt1 = t1.squaredNorm() ;

    if( isSmall( nt1 ) )
    {
        t1 = Vec3x ( 0, -normal[2], normal[1] ) ;
        t1.normalize() ;
    } else {
        t1 /= std::sqrt( nt1 ) ;
    }

    E.col(1) = t1 ;
    E.col(2) = normal.cross( t1 ).normalized() ;
}

void ProximityCollision::swapIfNecessary()
{
    if ( objects.first.globalIndex == -1 || (
            objects.second.globalIndex != -1
            && objects.second.globalIndex < objects.first.globalIndex ) )
    {
        normal = -normal;
        force = -force;
        std::swap( objects.first, objects.second );
    }
}

void ProximityCollision::print( std::ostream& os ) const
{
    os.precision( std::numeric_limits< double >::digits10 + 1);
    Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");

    os << "[PC:] S1: " << objects.first.globalIndex << " e1: " << objects.first.vertex << " a1: " << objects.first.abscissa
            << " S2: " << objects.second.globalIndex << " e1: " << objects.second.vertex << " a1: " << objects.second.abscissa
            << " N: " << normal.transpose().format(HeavyFmt) << " Rv: " << objects.first.freeVel.transpose().format(HeavyFmt) << "\n";
}

void ProximityCollision::printShort( std::ostream& os ) const
{
    os << "[" << objects.first.globalIndex << "." << objects.first.vertex << "|" <<
        objects.second.globalIndex << "." << objects.second.vertex << "] A: " << objects.first.abscissa << ", " << objects.second.abscissa
        << " N: " << normal.transpose() << " Rv: " << objects.first.freeVel.transpose();
}

bool ProximityCollision::operator< ( const ProximityCollision &rhs ) const
{
    if( objects.first.globalIndex == rhs.objects.first.globalIndex ){
        if ( objects.second.globalIndex == rhs.objects.second.globalIndex ){
            if( objects.first.vertex == rhs.objects.first.vertex ){
                return objects.second.vertex < rhs.objects.second.vertex;
            }
            else{
                return objects.first.vertex < rhs.objects.first.vertex;
            }
        }
        else{
            return objects.second.globalIndex < rhs.objects.second.globalIndex;
        }
    }
    else{
        return objects.first.globalIndex < rhs.objects.first.globalIndex;
    }

}