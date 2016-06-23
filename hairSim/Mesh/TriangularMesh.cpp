#include "TriangularMesh.h"

TriMesh::TriMesh( const MeshScriptingController* controller ) :
        m_associatedController( controller )
{

}

void TriMesh::tagBoundaryEdge( unsigned f, int firstVertex, int secondVertex )
{
    TriangularFace& face = m_faces[ f ] ;

    for( short k = 0 ; k < 3 ; ++ k )
    {
        const short next = ( k + 1 ) % 3 ;
        const short prev = ( k + 2 ) % 3 ;
        if( firstVertex == face.idx[k] )
        {
            if( secondVertex == face.idx[ next ] )
            {
                face.boundaryEdges |= ( 1<<k ) ;
            } else if ( secondVertex == face.idx[ prev ]  )
            {
                face.boundaryEdges |= ( 1<<prev ) ;
            }
        }
    }
}

void TriMesh::clear()
{
    m_vertices.clear();
    m_displacements.clear();
    m_faces.clear();
}

TriMesh::~TriMesh()
{
    // TODO Auto-generated destructor stub
}
