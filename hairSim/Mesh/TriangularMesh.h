#ifndef TRIANGULARMESH_HH_
#define TRIANGULARMESH_HH_

#include "../Core/Definitions.hh"
class TriMeshController;

struct TriangularFace
{ 
    int idx[3];
    int boundaryEdges ;

    TriangularFace( int un, int deux, int trois )
        : boundaryEdges( 0 )
    {
        idx[0] = un;
        idx[1] = deux;
        idx[2] = trois;
    }

    bool edgeOnBoundary( const short edge ) const
    {
        return ( boundaryEdges & (1 << edge) ) ;
    }

    bool hasBoundaryEdges() const
    {
        return boundaryEdges ;
    }
};

class TriMesh
{
public:
    TriMesh( const TriMeshController* controller = NULL );
    virtual ~TriMesh();

    void addVertex( const Vec3x& point )
    {
        m_vertices.push_back( point );
        m_displacements.push_back( Vec3x() );
    }

    Vec3x getVertex( size_t vtx ) const
    {
        return m_vertices[vtx];
    }

    void setVertex( size_t vtx, const Vec3x& point )
    {
        m_vertices[vtx] = point;
    }

    Vec3x getDisplacement( size_t vtx ) const
    {
        return m_displacements[vtx];
    }

    void setDisplacement( size_t vtx, const Vec3x& disp )
    {
        m_displacements[vtx] = disp;
    }

    const TriangularFace& getFace( unsigned f ) const
    {
        return m_faces[f];
    }

    void addFace( int un, int deux, int trois )
    {
        m_faces.push_back( TriangularFace( un, deux, trois ) );
    }

    void tagBoundaryEdge( unsigned f, int firstVertex, int secondVertex ) ;

    const std::vector<TriangularFace>& getFaces() const
    {
        return m_faces;
    }

    const TriMeshController* controller() const
    {
        return m_controller;
    }

    void setController( const TriMeshController* controller )
    {
        m_controller = controller;
    }

    size_t nv() const
    {
        return m_vertices.size();
    }

    size_t nf() const
    {
        return m_faces.size();
    }

    Vec3xArray &vertices()
    {
        return m_vertices;
    }

    void clear();

private:
    Vec3xArray m_vertices;
    Vec3xArray m_displacements;
    std::vector<TriangularFace> m_faces;
    const TriMeshController* m_controller;
};

#endif /* TRIANGULARMESH_HH_ */
