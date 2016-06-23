#ifndef STRAND_RENDERER_H
#define STRAND_RENDERER_H

#include "RenderBase.h"
#include "Color.h"
#include "OpenGLDecl.h"
#include "../Utils/ThreadUtils.h"
#include <set>

class ElementProxy;
class TwistEdgeHandler;
class ElasticStrand;

class StrandRenderer
{
public:

    enum DrawMode
    {
        LINES, QUADS, QUADS_SHADED, QUADS_BLENDED
    };

    enum DrawForce
    {
        NONE,
        GRAVITATION,
        BENDING,
        BENDING_VISCOUS,
        STRETCHING,
        STRETCHING_VISCOUS,
        TWISTING,
        TWISTING_VISCOUS
    };

    struct QuadData
    {
        std::vector<GLfloat> m_quadVertices;
        std::vector<GLfloat> m_quadColors;
        std::vector<GLfloat> m_quadNormals;
        std::vector<GLuint> m_quadIndices;
    };

    StrandRenderer();

    void render( ElasticStrand* strand, const int& w, const int& h, const int& label, const bool& ct );
    void addVert( Vec3x& vert );
    void addArrow( Vec3x& start, Vec3x& end );
    static Vec3x calculateObjectCenter( ElasticStrand* strand );
    static Scalar calculateObjectBoundingRadius( ElasticStrand* strand, const Vec3x& center );

    const Eigen::Matrix4f &transformationMatrix() const
    {
        return m_transformationMatrix ;
    }

    Eigen::Matrix4f &transformationMatrix()
    {
        return m_transformationMatrix ;
    }

    void pushTransformationMatrix() const ;
    void pushInverseTransformationMatrix() const ;
    void popTransformationMatrix() const ;

/////////
    void labelScene( int windowWidth, int windowHeight, int label );
    void drawProxies( std::vector<ElementProxy* >* proxies, TwistEdgeHandler* teh );
    void renderProxies();

    VecXx displaceVec;
    std::vector<Vec3x> verts;
    std::vector< std::pair<Vec3x, Vec3x> > arrows;
    std::vector< std::pair<Vec3x, Vec3x> > startProxies;
    std::vector< std::pair<Vec3x, Vec3x> > endProxies;
    std::vector< std::pair<Vec3x, int> > vertLabels;
    std::vector< std::pair<Vec3x, double> > twistLabels;
    int m_wWidth;
    int m_wHeight;
    int m_label;
////////

protected:

    void drawSmoothStrand();
    void drawVertices( int flag = GL_LINE_STRIP ) const;
    void drawContacts() const;
    void drawArrows() const;

    template<typename ForceT>
    void drawForce(); 
    void drawForceVec( VecXx& f );
    void computeQuads( QuadData &quads ) const ;
    static void drawQuads( const QuadData &quads ) ;

    ElasticStrand* m_strand;
    std::map<const char*, Color> m_palette;
    Eigen::Matrix4f m_transformationMatrix ;
    double m_strandRadius;
    DrawMode m_drawMode;
    DrawForce m_drawForce;
    QuadData m_smoothStrand ;

};

#endif // RODRENDERER_HH
