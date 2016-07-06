#include "StrandRenderer.h"

#include "../Strand/ElasticStrand.h"
#include "../Strand/StrandState.h"
#include "../Strand/ElasticStrandUtils.h"
#include "../Forces/GravitationForce.hh"
#include "../Forces/BendingForce.hh"
#include "../Forces/StretchingForce.hh"
#include "../Forces/TwistingForce.hh"
#include "../Forces/ForceAccumulator.hh"

#include "../Collision/ElementProxy.h"
#include "../Utils/StringUtils.h"
#include "RenderUtils.h"

#define RADPERDEG 0.0174533
#define NUM_PER_LOCK 31

StrandRenderer::StrandRenderer():
        m_strand( NULL ), //
        m_drawMode( SHADED ), // LINES, SHADED, BLENDED
        m_drawForce( NONE ) // NONE, GRAVITATION, BENDING, STRETCHING, TWISTING
{
    m_transformationMatrix.setIdentity();
    m_palette["random"] = Color(( (float) (rand() % 255) ) / 255., ( (float) (rand() % 255) ) / 255., ( (float) (rand() % 255) ) / 255.);
    m_palette["orange"] = Color( 189, 73, 50 );
    m_palette["yellow"] = Color( 255, 211, 78 );
    m_palette["green"] = Color( 16, 91, 99 );
    m_palette["force"] = Color( 0, 0, 255 );
}

void StrandRenderer::render( ElasticStrand* strand, const int& w, const int& h, const int& label, const bool& ct )
{
    m_strand = strand;
    m_wWidth = w;
    m_wHeight = h;
    m_label = (label + 2) % 3;

    if( ct ) m_strandRadius = m_strand->collisionParameters().collisionRadius( 0 );
    else m_strandRadius = m_strand->m_physicsRadius;

    glPushAttrib( GL_COLOR_BUFFER_BIT );
    pushTransformationMatrix();

    drawSmoothStrand();
    drawVertices();
    drawContacts();
    drawArrows();
    // renderProxies();

    switch ( m_drawForce )
    {
    case GRAVITATION:
        drawForce< GravitationForce >();
        break;
    case BENDING:
        drawForce< BendingForce<> >();
        break;
    case BENDING_VISCOUS:
        drawForce< BendingForce<Viscous> >();
        break;
    case STRETCHING:
        drawForce< StretchingForce<> >();
        break;
    case STRETCHING_VISCOUS:
        drawForce< StretchingForce<Viscous> >();
        break;
    case TWISTING:
        drawForce< TwistingForce<> >();
        break;
    case TWISTING_VISCOUS:
        drawForce< TwistingForce<Viscous> >();
        break;
    case NONE:
    default:
        break;
    }

    if( displaceVec.size() > 1 && !displaceVec.isZero() ){
        drawForceVec( displaceVec );
    }

    popTransformationMatrix();
    glPopAttrib();
}

void StrandRenderer::addVert( Vec3& vert )
{
    verts.push_back( vert );
}

void StrandRenderer::addArrow( Vec3& start, Vec3& end ){
    arrows.push_back( std::pair< Vec3, Vec3>( start, end ) );
}

void StrandRenderer::drawContacts() const
{
    glPushAttrib( GL_ALL_ATTRIB_BITS );
    glPointSize( 10 );
    glColor4f( 120.0 / 255 , 0.0 , 1.0, 0.5 ); // purple
    glBegin( GL_POINTS );
    for( unsigned v = 0; v < verts.size(); ++v ){
        glVertex3dv( verts[v].data() );
    }
    glEnd();
    glPopAttrib();
}

void StrandRenderer::drawVertices() const
{
    pushMode();

    int numperlock = NUM_PER_LOCK;
    const char* c = NULL;
    if( m_strand->getGlobalIndex() / numperlock == 0 ) c = "orange";
    if( m_strand->getGlobalIndex() / numperlock == 1 ) c = "yellow";
    if( m_strand->getGlobalIndex() / numperlock > 1 ) c = "green";
    const Color &co = m_palette.find( c )->second;
    glColor3dv( co.data() );

    for( unsigned v = 1; v < m_strand->getNumVertices() - 1; ++v )
    {
        renderSphere( m_strand->getVertex( v ), m_strandRadius );
    }

    popMode();
}

void Arrow( GLdouble x1, GLdouble y1, GLdouble z1, GLdouble x2, GLdouble y2, GLdouble z2, GLdouble D )
{
    double x=x2-x1;
    double y=y2-y1;
    double z=z2-z1;
    double L=sqrt(x*x+y*y+z*z);

    GLUquadricObj *quadObj;

    glPushMatrix ();
        glColor3f( 0., 0., 1. );

        glTranslated(x1,y1,z1);

        if((x!=0.)||(y!=0.)) {
            glRotated(atan2(y,x)/RADPERDEG,0.,0.,1.);
            glRotated(atan2(sqrt(x*x+y*y),z)/RADPERDEG,0.,1.,0.);
        } else if (z<0){
            glRotated(180,1.,0.,0.);
        }

        glTranslatef(0,0,L-4*D);

        quadObj = gluNewQuadric ();
        gluQuadricDrawStyle (quadObj, GLU_FILL);
        gluQuadricNormals (quadObj, GLU_SMOOTH);
        gluCylinder(quadObj, 2*D, 0.0, 4*D, 32, 1);
        gluDeleteQuadric(quadObj);

        quadObj = gluNewQuadric ();
        gluQuadricDrawStyle (quadObj, GLU_FILL);
        gluQuadricNormals (quadObj, GLU_SMOOTH);
        gluDisk(quadObj, 0.0, 2*D, 32, 1);
        gluDeleteQuadric(quadObj);

        glTranslatef(0,0,-L+4*D);

        quadObj = gluNewQuadric ();
        gluQuadricDrawStyle (quadObj, GLU_FILL);
        gluQuadricNormals (quadObj, GLU_SMOOTH);
        gluCylinder(quadObj, D, D, L-4*D, 32, 1);
        gluDeleteQuadric(quadObj);

        quadObj = gluNewQuadric ();
        gluQuadricDrawStyle (quadObj, GLU_FILL);
        gluQuadricNormals (quadObj, GLU_SMOOTH);
        gluDisk(quadObj, 0.0, D, 32, 1);
        gluDeleteQuadric(quadObj);

    glPopMatrix ();
}

void StrandRenderer::drawArrows() const
{
    for (int i = 0; i < arrows.size(); ++i)
    {
        Vec3 start = arrows[i].first;
        Vec3 end = arrows[i].second;
        Arrow( start[0], start[1], start[2], end[0], end[1], end[2], 0.02 );
    }
}

void StrandRenderer::drawCylinders()
{
    int numperlock = NUM_PER_LOCK;
    const char* c = NULL ;
    if( m_strand->getGlobalIndex() / numperlock == 0 ) c = "orange";
    if( m_strand->getGlobalIndex() / numperlock == 1 ) c = "yellow";
    if( m_strand->getGlobalIndex() / numperlock > 1 ) c = "green";
    const Color &co = m_palette.find( c )->second;
    glColor3dv( co.data() );

    for( int v = 0; v < m_strand->getNumVertices() - 1; ++v )
    {
        const Vec3& a = m_strand->getVertex( v );
        const Vec3& b = m_strand->getVertex( v + 1 );
        renderCylinder( a, b, m_strandRadius );
    }
}

void StrandRenderer::drawSmoothStrand()
{
    pushMode();
    drawCylinders();
    popMode();
}

template<typename ForceT>
void StrandRenderer::drawForce()
{
    static const double SCALE = 2.0;

    glPushAttrib( GL_COLOR );

    VecXx force( m_strand->getNumVertices() * 4 - 1 );

    // Here is a clever/ugly hack: simply doing accumulateCurrentF would give zero for viscous forces
    // because they take the current state as "rest shape".
    // Actually, if you try to call accumulateCurrentF on a dissipative force it won't compile, see why in ForceAccumulator
    const_cast<ElasticStrand&>( *m_strand ).swapStates();
    ForceAccumulator<ForceT>::accumulateFuture( force, const_cast<ElasticStrand&>( *m_strand) );
    const_cast<ElasticStrand&>( *m_strand ).swapStates();

    glLineWidth( 2 );
    glBegin( GL_LINES );
    glColor3dv( m_palette["force"].data() );

    for ( int vtx = 0; vtx < m_strand->getNumVertices(); ++vtx )
    {
        const Eigen::Vector3d x0 = m_strand->getVertex( vtx );
        glVertex3dv( x0.data() );
        const Eigen::Vector3d x1 = x0 + SCALE * force.segment<3>( 4 * vtx );
        glVertex3dv( x1.data() );
    }

    glEnd();
    glPopAttrib();
}

void StrandRenderer::drawForceVec( VecXx& f )
{
    static const double SCALE = 1.0;
    glPushAttrib( GL_COLOR );

    glLineWidth( 2 );
    glBegin( GL_LINES );
    glColor3dv( m_palette["force"].data() );

    for ( int vtx = 0; vtx < m_strand->getNumVertices(); ++vtx )
    {
        const Eigen::Vector3d x0 = m_strand->getVertex( vtx );
        glVertex3dv( x0.data() );
        const Eigen::Vector3d x1 = x0 + SCALE * f.segment<3>( 4 * vtx );
        glVertex3dv( x1.data() );
    }

    glEnd();
    glPopAttrib();
}

void renderBitmapString( float x, float y, float z, void *font, std::string s ) 
{
    glRasterPos3f(x, y, z);
    for( std::string::iterator i = s.begin(); i != s.end(); ++i )
    {
        char c = *i;
        glutBitmapCharacter(font, c);
    }
}

void StrandRenderer::labelScene( int windowWidth, int windowHeight, int label )
{
    GLint viewport[4];
    GLdouble modelview[16];
    GLdouble viewVector[3];
    GLdouble projection[16];

    GLdouble winX, winY, winZ; // 2D point
    GLdouble posX, posY, posZ; // 3D point
    posX = 0.0;
    posY = 0.0;
    posZ = -1.0; // the display is the same if posZ=1 which should not be the case

    //get the matrices
    glGetDoublev( GL_MODELVIEW_MATRIX, modelview );
    viewVector[0]=modelview[8];
    viewVector[1]=modelview[9];
    viewVector[2]=modelview[10];
    glGetDoublev( GL_PROJECTION_MATRIX, projection );
    glGetIntegerv( GL_VIEWPORT, viewport );

    Vec3 x1;
    
    if( label < 2 )
    {

    for( unsigned i = 0; i < vertLabels.size(); ++i )
    {
        x1 = vertLabels[i].first;
        posX = x1[0];
        posY = x1[1];
        posZ = x1[2];
        
        int res = gluProject(posX,posY,posZ,modelview,projection,viewport,&winX,&winY,&winZ);
        
        // if( res && viewVector[0]*posX + viewVector[1]*posY + viewVector[2]*posZ > 0 ){
            
            glMatrixMode(GL_PROJECTION);
            glPushMatrix();
            glLoadIdentity();
            gluOrtho2D(0, windowWidth, 0, windowHeight);
            glMatrixMode(GL_MODELVIEW);
            glPushMatrix();
            glLoadIdentity();
            // assert( renderingutils::checkGLErrors() );
            
            glColor3f( 1.0, 1.0, 1.0 );

            if( (vertLabels[i].second) < 0 ){
                renderBitmapString( winX, winY, 0.0, GLUT_BITMAP_HELVETICA_12, "v" + toString( -1* (vertLabels[i].second) ) );                
            }
            else{
                renderBitmapString( winX, winY, 0.0, GLUT_BITMAP_HELVETICA_12, "E" + toString( vertLabels[i].second ) );                                
            }

            glMatrixMode(GL_PROJECTION);
            glPopMatrix();
            glMatrixMode(GL_MODELVIEW);
            glPopMatrix();
        // }
        
    } 
    }
    
    if( label == 0 ){
        for( unsigned i = 0; i < twistLabels.size(); i += 3 )
        {
            x1 = twistLabels[i].first;
            posX = x1[0];
            posY = x1[1];
            posZ = x1[2];
            
            int res = gluProject(posX,posY,posZ,modelview,projection,viewport,&winX,&winY,&winZ);
            
            // if( res && viewVector[0]*posX + viewVector[1]*posY + viewVector[2]*posZ < 0 ){
                
                glMatrixMode(GL_PROJECTION);
                glPushMatrix();
                glLoadIdentity();
                gluOrtho2D(0, windowWidth, 0, windowHeight);
                glMatrixMode(GL_MODELVIEW);
                glPushMatrix();
                glLoadIdentity();
                // assert( renderingutils::checkGLErrors() );

                glColor3f( 0.0, 0.0, 0.0 );

                renderBitmapString( winX, winY, 0.0, GLUT_BITMAP_HELVETICA_12, "         coplanar:" + toString( twistLabels[i].second  ) );
                renderBitmapString( winX, winY - 12, 0.0, GLUT_BITMAP_HELVETICA_12, "         intersection:" + toString( twistLabels[i + 1].second ) );
                renderBitmapString( winX, winY - 24, 0.0, GLUT_BITMAP_HELVETICA_12, "         twistAngle:" + toString( twistLabels[i + 2].second * (180.0/3.14159) ) );

                glMatrixMode(GL_PROJECTION);
                glPopMatrix();
                glMatrixMode(GL_MODELVIEW);
                glPopMatrix();
            // }
            
        } 
    }      
        
}

void StrandRenderer::renderProxies( )
{
    glPushAttrib( GL_COLOR );

    glLineWidth( 2 );
    glBegin( GL_LINES );
    glColor3dv( m_palette["velocity"].data() );
    for( unsigned a = 0; a < endProxies.size(); ++a )
    {
            glVertex3dv( endProxies[a].first.data() );
            glVertex3dv( endProxies[a].second.data() );
    }

    glEnd();
    glPopAttrib();

    labelScene( m_wWidth, m_wHeight, m_label );
}

void StrandRenderer::drawProxies( std::vector<ElementProxy* >* proxies, TwistEdgeHandler* teh )
{
    for( unsigned e = 0; e < proxies->size(); ++e )
    {
        TwistEdge* edge = dynamic_cast<TwistEdge*>( (*proxies)[e] );
        if( edge && edge->intersectionTwists() > 0 ){
            Vec3 x0, x1;
            teh->getEdgeVerts( edge, false, x0, x1 );
            endProxies.push_back( std::pair< Vec3, Vec3>( x0, x1 ) );

            Vec3 midpoint = 0.5 *( x0 + x1 );
            vertLabels.push_back( std::pair< Vec3, int >( midpoint, edge->uniqueID ) );

            if( edge->isTwistedBand ){
                twistLabels.push_back( std::pair< Vec3, double >(midpoint, edge->coplanarTwists() ) );
                twistLabels.push_back( std::pair< Vec3, double >(midpoint, edge->intersectionTwists() ) );
                twistLabels.push_back( std::pair< Vec3, double >(midpoint, edge->currAngle() ) );
            }

            // if( edge->endpoints.second->incidentEdges.second == NULL ) vertLabels.push_back( std::pair< Vec3, int >(x0, -1 * (edge->endpoints.second->uniqueId) ) );
        }
    }
}

Vec3 StrandRenderer::calculateObjectCenter( ElasticStrand* strand )
{
    Vec3 center = Vec3::Zero();

    const unsigned numVertices = strand->getNumVertices();
    for( unsigned vtx = 0; vtx < numVertices; ++vtx )
    {
        const Vec3& vertex = strand->getVertex( vtx );
        center += vertex;
    }

    center /= numVertices;
    return center;
}

Scalar StrandRenderer::calculateObjectBoundingRadius( ElasticStrand* strand, const Vec3& center )
{
    Scalar radius = 0.0;

    const unsigned numVertices = strand->getNumVertices();
    for ( unsigned vtx = 0; vtx < numVertices; ++vtx )
    {
        const Vec3& vertex = strand->getVertex( vtx );
        radius = std::max( radius, ( vertex - center ).norm() );
    }
    return radius;
}

void StrandRenderer::pushMode() const
{
    if ( m_drawMode >= SHADED )
    {
        glEnable( GL_LIGHTING );
        glEnable( GL_COLOR_MATERIAL );
        glColorMaterial( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE );

        // Hack to make it work with maya's default two-sided lighting
        glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE );

        if ( m_drawMode == BLENDED )
        {
            glEnable( GL_BLEND );
            glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
        }
    }
}

void StrandRenderer::popMode() const
{
    if ( m_drawMode >= SHADED )
    {
        glDisable( GL_COLOR_MATERIAL );
        glDisable( GL_LIGHTING );

        if ( m_drawMode == BLENDED )
        {
            glDisable( GL_BLEND );
        }
    }
}

void StrandRenderer::pushTransformationMatrix() const
{
    glPushMatrix();
    glMultMatrixf( m_transformationMatrix.data() );
}

void StrandRenderer::pushInverseTransformationMatrix() const
{
    glPushMatrix();
    Eigen::Matrix4f inv( m_transformationMatrix.inverse() );
    glMultMatrixf( inv.data() );
}

void StrandRenderer::popTransformationMatrix() const
{
    glPopMatrix();
}
