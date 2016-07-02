#include "ElementProxy.h"
#include "CollisionUtils/BVH.hh"
#include "../Utils/StringUtils.h"

#define SQRT_2 1.41421356237

std::ostream& operator<<( std::ostream& os, const ElementProxy& elem )
{
    elem.print( os );
    return os;
}

void EdgeProxy::print( std::ostream& os ) const
{
    os << "edge: " << &m_strand << ' ' << m_vertexIndex << ": "
        << m_strand.getFutureVertex( m_vertexIndex ).format( EIGEN_VECTOR_IO_FORMAT ) << " --- "
        << m_strand.getFutureVertex( m_vertexIndex + 1 ).format( EIGEN_VECTOR_IO_FORMAT );
}

void FaceProxy::print( std::ostream& os ) const
{
    const TriMesh* const mesh = getMesh();

    os << "triangle: " << mesh << '\n';

    os << "past position:\n";
    os
            << m_face.idx[0]
            << ": "
            << ( mesh->getVertex( m_face.idx[0] ) - mesh->getDisplacement( m_face.idx[0] ) ).format(
                    EIGEN_VECTOR_IO_FORMAT ) << '\n';
    os
            << m_face.idx[1]
            << ": "
            << ( mesh->getVertex( m_face.idx[1] ) - mesh->getDisplacement( m_face.idx[1] ) ).format(
                    EIGEN_VECTOR_IO_FORMAT ) << '\n';
    os
            << m_face.idx[2]
            << ": "
            << ( mesh->getVertex( m_face.idx[2] ) - mesh->getDisplacement( m_face.idx[2] ) ).format(
                    EIGEN_VECTOR_IO_FORMAT ) << '\n';

    os << "current position:\n";
    os << m_face.idx[0] << ": " << mesh->getVertex( m_face.idx[0] ).format( EIGEN_VECTOR_IO_FORMAT )
            << '\n';
    os << m_face.idx[1] << ": " << mesh->getVertex( m_face.idx[1] ).format( EIGEN_VECTOR_IO_FORMAT )
            << '\n';
    os << m_face.idx[2] << ": "
            << mesh->getVertex( m_face.idx[2] ).format( EIGEN_VECTOR_IO_FORMAT );
}

void TwistEdge::computeBoundingBox( BBoxType& boundingBox, bool statique, TwistEdgeHandler* teh ) const
{
    if( !isTwistedBand ){
        CylinderProxy::computeBoundingBox( boundingBox, statique );
        return;
    }

    boundingBox.reset();
    Vec3f min, max;
    Vec3 vA, vB;
    teh->getEdgeVerts( this, true, vA, vB );

    boundingBox.insert( vA.cast< float >() );
    boundingBox.insert( vB.cast< float >() );

    if ( !statique ) 
    {        
        teh->getEdgeVerts( this, false, vA, vB );
        boundingBox.insert( vA.cast< float >() );
        boundingBox.insert( vB.cast< float >() );
    }

    static const Vec3f unit = Vec3f::Ones();
    min = boundingBox.min - m_radius * unit;
    max = boundingBox.max + m_radius * unit;

    boundingBox.insert( min );
    boundingBox.insert( max );
}

int edgeCounter = 0;
TwistEdge::TwistEdge( ElasticStrand& strand, int vertexIndex, ImplicitStepper* stepper ):
    CylinderProxy( strand, vertexIndex, stepper ),
    isTwistedBand( false ), 
    flagged( false ),
    traversed( false ),
    prev( NULL ),
    next( NULL )
{
    uniqueID = edgeCounter++;
    parents.first = NULL;
    parents.second = NULL;
}

TwistEdge::TwistEdge( TwistEdge* first, TwistEdge* second ):
    CylinderProxy( first->m_strand, first->m_vertexIndex, first->m_implicit ),
    isTwistedBand( true ),
    flagged( false ),
    traversed( false ),
    prev( NULL ),
    next( NULL )
{
    TwistEdge* alpha = first;
    TwistEdge* beta = second;

    uniqueID = edgeCounter++;

    // std::cout << "Creating tunneling from first: " << first->uniqueID << " second: " << second->uniqueID << " with edgeCounter: " << edgeCounter << " my uniqueID: " << uniqueID << std::endl;

    parents.first = alpha;
    parents.second = beta;

    // intersections.push( new TwistIntersection( angleAtTimeOfCollision ) );

    alpha->children.push_back( std::pair< int, TwistEdge* >( beta->uniqueID, this ) );
    beta->children.push_back( std::pair< int, TwistEdge* >( alpha->uniqueID, this ) );
}

TwistEdge::TwistEdge( TwistEdge* first, TwistEdge* second, int& uID ):
    CylinderProxy( first->m_strand, first->m_vertexIndex, first->m_implicit ),
    isTwistedBand( true ),
    flagged( false ),
    traversed( false ),
    prev( NULL ),
    next( NULL )
{
    TwistEdge* alpha = first;
    TwistEdge* beta = second;

    uniqueID = uID;

    // std::cout << "reinstating TwistEdge: " << uniqueID << " between " << first->uniqueID << " second: " << second->uniqueID << std::endl;

    parents.first = alpha;
    parents.second = beta;

    // intersections.push( new TwistIntersection( angleAtTimeOfCollision ) );

    alpha->children.push_back( std::pair< int, TwistEdge* >( beta->uniqueID, this ) );
    beta->children.push_back( std::pair< int, TwistEdge* >( alpha->uniqueID, this ) );
}

void TwistEdge::printIntersections( std::stack< TwistIntersection* >& intersections, std::ofstream& out )
{
    if( intersections.empty() )
    {
        // out << std::endl;
        return;
    }
    TwistIntersection* x = intersections.top();
    intersections.pop();
    printIntersections( intersections, out );
    intersections.push( x );
    serializeVarHex( x->originalTwistAngle, out );
    serializeVarHex( x->currAngle, out );
    serializeVarHex( x->coplanarTwists, out );

    // out << x->originalTwistAngle << " " << x->currAngle << " " << x->coplanarTwists << " ";
}


Scalar TwistEdge::currAngle() const
{
    if( intersections.size() == 0 ) return 0;
    return intersections.top()->currAngle;
}

int TwistEdge::coplanarTwists() const
{
    if( intersections.size() == 0 ) return 0;
    return intersections.top()->coplanarTwists;
}

int TwistEdge::intersectionTwists() const
{
    return intersections.size();
}

Scalar TwistEdge::originalTwistAngle() const
{
    if( intersections.size() == 0 ) return 0;
    return intersections.top()->originalTwistAngle;
}


