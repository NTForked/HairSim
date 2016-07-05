#ifndef STRAND_RENDERER_H
#define STRAND_RENDERER_H

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
        LINES, SHADED, BLENDED
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

    StrandRenderer();

    void render( ElasticStrand* strand, const int& w, const int& h, const int& label, const bool& ct );
    void addVert( Vec3& vert );
    void addArrow( Vec3& start, Vec3& end );
    static Vec3 calculateObjectCenter( ElasticStrand* strand );
    static Scalar calculateObjectBoundingRadius( ElasticStrand* strand, const Vec3& center );

    const Eigen::Matrix4f &transformationMatrix() const
    { return m_transformationMatrix; }

    Eigen::Matrix4f &transformationMatrix()
    { return m_transformationMatrix; }

    void pushTransformationMatrix() const;
    void pushInverseTransformationMatrix() const;
    void popTransformationMatrix() const;

/////////
    void labelScene( int windowWidth, int windowHeight, int label );
    void drawProxies( std::vector<ElementProxy* >* proxies, TwistEdgeHandler* teh );
    void renderProxies();

    VecXx displaceVec;
    std::vector<Vec3> verts;
    std::vector< std::pair<Vec3, Vec3> > arrows;
    std::vector< std::pair<Vec3, Vec3> > startProxies;
    std::vector< std::pair<Vec3, Vec3> > endProxies;
    std::vector< std::pair<Vec3, int> > vertLabels;
    std::vector< std::pair<Vec3, double> > twistLabels;
    int m_wWidth;
    int m_wHeight;
    int m_label;
////////

protected:

    void drawSmoothStrand();
    void drawVertices() const;
    void drawContacts() const;
    void drawArrows() const;

    template<typename ForceT>
    void drawForce(); 
    void drawForceVec( VecXx& f );
    void drawCylinders();

    void pushMode() const;
    void popMode() const;

    ElasticStrand* m_strand;
    std::map<const char*, Color> m_palette;
    Eigen::Matrix4f m_transformationMatrix ;
    double m_strandRadius;
    DrawMode m_drawMode;
    DrawForce m_drawForce;

};

#endif // RODRENDERER_HH
