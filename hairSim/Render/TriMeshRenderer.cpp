#include "TriangleMeshRenderer.hh"

#include "../Render/Color.hh"
#include "../Render/OpenGLDecl.hh"

TriangleMeshRenderer::TriangleMeshRenderer( const TriangularMesh& mesh )
: m_mesh(mesh)
, m_mode(FLAT)
{}

void TriangleMeshRenderer::render()
{
    if( m_mode == FLAT )
    {
        glEnable(GL_LIGHTING);
        
        GLfloat gray[] = {0.8,0.8,0.8,1.0};
        glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,gray);
        
        // Render all faces
        glBegin(GL_TRIANGLES);
        
        for( auto fit = m_mesh.getFaces().begin(); fit != m_mesh.getFaces().end(); ++fit )
        {
            std::vector<Vec3d> v;
            for( int i=0; i<3 ;++i ){
                v.push_back(m_mesh.getVertex((*fit).idx[i]));
            }
            
            // Compute a normal for the face
            Vec3d e1 = v[1]-v[0];
            Vec3d e2 = v[2]-v[0];
            Vec3d n = e1.cross(e2);
            if( n.norm() != 0 ) n.normalize();
            glNormal3f(n.x(),n.y(),n.z());
            glVertex3f(v[0].x(),v[0].y(),v[0].z());
            glVertex3f(v[1].x(),v[1].y(),v[1].z());
            glVertex3f(v[2].x(),v[2].y(),v[2].z());
        }
        glEnd();
        
        glDisable(GL_LIGHTING);
    }
    else if( m_mode == DBG )
    {
        glDisable(GL_LIGHTING);
                
        // Render all faces
        glBegin(GL_TRIANGLES);
        OpenGL::color(Color(255,0,0));

        for( auto fit = m_mesh.getFaces().begin(); fit != m_mesh.getFaces().end(); ++fit )
        {
            for( int i=0; i<3 ; ++i )
            {
                OpenGL::vertex( m_mesh.getVertex((*fit).idx[i]) );
            }
        }
        glEnd();
        
        // Render all vertices
        glPointSize(5);
        glBegin(GL_POINTS);
        OpenGL::color(Color(0,0,0));
        for( unsigned i = 0; i < m_mesh.nv(); ++i )
        {
            OpenGL::vertex(m_mesh.getVertex(i));
        }
        
        glEnd();
        
        glEnable(GL_LIGHTING);
    }
}

Vec3d TriangleMeshRenderer::calculateObjectCenter()
{
    Vec3d center( 0.0,0.0,0.0 );
    for( unsigned i = 0; i < m_mesh.nv(); ++i )
    {
        center += m_mesh.getVertex(i);
    }
    
    if( m_mesh.nv() != 0 ) center /= ((double)m_mesh.nv());
    
    return center;
}

double TriangleMeshRenderer::calculateObjectBoundingRadius( const Vec3d& center )
{
    Scalar radius = 0.0;
    for( unsigned i = 0; i < m_mesh.nv(); ++i )
    {
        radius = std::max(radius, (m_mesh.getVertex(i)- center).norm());
    }
    
    return radius;
}

