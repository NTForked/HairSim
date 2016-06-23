#include "ElementProxy.hh"
#include "BVH.hh"
#include "../Dynamic/Config.hh"

#include "../Dynamic/MeshScriptingController.hh"

#define SQRT_2 1.41421356237

std::ostream& operator<<( std::ostream& os, const ElementProxy& elem )
{
    elem.print( os );
    return os;
}

void EdgeProxy::print( std::ostream& os ) const
{
    os << "edge: " << &m_strand << ' ' << m_vertexIndex << ": "
        << m_strand.getVertex( m_vertexIndex ).format( EIGEN_VECTOR_IO_FORMAT ) << " --- "
        << m_strand.getVertex( m_vertexIndex + 1 ).format( EIGEN_VECTOR_IO_FORMAT );
}

void FaceProxy::print( std::ostream& os ) const
{
    const TriangularMesh* const mesh = getMesh();

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

    std::cerr << "this stuff below is definitely wrong, do not need ortho" << std::endl;

    boundingBox.reset();
    Vec3x min, max;
    Vec3x vA, vB;
    teh->getEdgeVerts( this, true, vA, vB );
    min = vA.cwiseMin(vB);
    max = vA.cwiseMax(vB);

    static const Vec3x unit = Vec3x::Ones();
    Vec3x ortho = (vB-vA).cross(unit); ortho.normalize();
    Vec3x ortho2 = unit.cross(vB-vA); ortho2.normalize();
    min += SQRT_2 * radius * ortho.cwiseMin(ortho2);
    max += SQRT_2 * radius * ortho.cwiseMax(ortho2);
    boundingBox.insert( min );
    boundingBox.insert( max );

    if ( !statique ) 
    {        
        teh->getEdgeVerts( this, false, vA, vB );
        min = vA.cwiseMin(vB);
        max = vA.cwiseMax(vB);
        ortho = (vB-vA).cross(unit); ortho.normalize();
        ortho2 = unit.cross(vB-vA); ortho2.normalize();
        min += SQRT_2 * radius * ortho.cwiseMin(ortho2);
        max += SQRT_2 * radius * ortho.cwiseMax(ortho2);
        boundingBox.insert( min );
        boundingBox.insert( max );
    }
}

int edgeCounter = 0;
TwistEdge::TwistEdge( ElasticStrand& strand, int vertexIndex, Scalar radius, ImplicitStepper* stepper ):
    EdgeProxy( strand, vertexIndex, stepper ),
    isTwistedBand( false ), 
    radius( radius ),
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
    EdgeProxy( first->m_strand, first->m_vertexIndex, first->m_implicit ),
    isTwistedBand( true ),
    flagged( false ),
    traversed( false ),
    prev( NULL ),
    next( NULL )
{
    TwistEdge* alpha = first;
    TwistEdge* beta = second;

    radius = 0.5 * (alpha->radius + beta->radius);
    uniqueID = edgeCounter++;

    // std::cout << "Creating tunneling from first: " << first->uniqueID << " second: " << second->uniqueID << " with edgeCounter: " << edgeCounter << " my uniqueID: " << uniqueID << std::endl;

    parents.first = alpha;
    parents.second = beta;

    // intersections.push( new TwistIntersection( angleAtTimeOfCollision ) );

    alpha->children.push_back( std::pair< int, TwistEdge* >( beta->uniqueID, this ) );
    beta->children.push_back( std::pair< int, TwistEdge* >( alpha->uniqueID, this ) );
}

TwistEdge::TwistEdge( TwistEdge* first, TwistEdge* second, int& uID ):
    EdgeProxy( first->m_strand, first->m_vertexIndex, first->m_implicit ),
    isTwistedBand( true ),
    flagged( false ),
    traversed( false ),
    prev( NULL ),
    next( NULL )
{
    TwistEdge* alpha = first;
    TwistEdge* beta = second;

    radius = 0.5 * (alpha->radius + beta->radius);
    uniqueID = uID;

    // std::cout << "reinstating TwistEdge: " << uniqueID << " between " << first->uniqueID << " second: " << second->uniqueID << std::endl;

    parents.first = alpha;
    parents.second = beta;

    // intersections.push( new TwistIntersection( angleAtTimeOfCollision ) );

    alpha->children.push_back( std::pair< int, TwistEdge* >( beta->uniqueID, this ) );
    beta->children.push_back( std::pair< int, TwistEdge* >( alpha->uniqueID, this ) );
}

template<typename T> void serializeVarHex( const T& var, std::ostream& output_stream )
{
    assert( output_stream.good() );
    T local_var = var;
    output_stream.precision( std::numeric_limits<double>::digits10 + 2);
    output_stream.flags( std::ios_base::fixed | std::ios_base::scientific );
    output_stream.write( reinterpret_cast<char*>( &local_var ), sizeof(T) );
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


