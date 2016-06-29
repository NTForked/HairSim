#ifndef TRI_MESH_RENDERER_H
#define TRI_MESH_RENDERER_H

#include "../Mesh/TriMesh.h"

/** Class that implements OpenGL rendering for triangle meshes. */
class TriMeshRenderer
{
public:
    
    enum DrawMode { DBG, FLAT, NONE };
    
    explicit TriMeshRenderer( const TriMesh& mesh );
    
    void render();
    
    DrawMode getMode() const { return m_mode; }
    void setMode(DrawMode mode) { m_mode = mode; }
    
    virtual Vec3d calculateObjectCenter();
    virtual Scalar calculateObjectBoundingRadius(const Vec3d& center);
    
protected:
    const TriMesh& m_mesh;
    DrawMode m_mode;
};


#endif 
