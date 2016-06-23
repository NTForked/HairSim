#ifndef TRI_MESH_RENDERER_H
#define TRI_MESH_RENDERER_H

#include "../Collision/TriMesh.h"

/** Class that implements OpenGL rendering for triangle meshes. */
class TriMeshRenderer
{
public:
    
    enum DrawMode { DBG, FLAT, NONE };
    
    explicit TriangleMeshRenderer( const TriangularMesh& mesh );
    
    void render();
    
    DrawMode getMode() const { return m_mode; }
    void setMode(DrawMode mode) { m_mode = mode; }
    
    virtual Vec3d calculateObjectCenter();
    virtual Scalar calculateObjectBoundingRadius(const Vec3d& center);
    
protected:
    const TriangularMesh& m_mesh;
    DrawMode m_mode;
};


#endif 
