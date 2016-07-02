#include "RenderUtils.h"

#define RADPERDEG 0.0174533

void renderSphere( const Vec3& p, GLdouble radius )
{
    GLUquadricObj *quadObj;
    glPushMatrix ();

        glTranslated( p[0], p[1], p[2] );

        quadObj = gluNewQuadric();
        gluQuadricDrawStyle( quadObj, GLU_FILL );
        gluSphere( quadObj, radius, 12, 12 );
        gluDeleteQuadric(quadObj);

    glPopMatrix ();
}

void renderCylinder( const Vec3& p, const Vec3& q, GLdouble radius )
{
    double x = q[0] - p[0];
    double y = q[1] - p[1];
    double z = q[2] - p[2];
    double L = sqrt( x*x + y*y + z*z );

    GLUquadricObj *quadObj;
    glPushMatrix ();

        glTranslated( p[0], p[1], p[2] );

        if( (x!=0.) || (y!=0.) )
        {
            glRotated( atan2( y, x ) / RADPERDEG, 0., 0., 1. );
            glRotated( atan2( sqrt( x*x + y*y ), z ) / RADPERDEG, 0., 1., 0.);
        } 
        else if( z < 0 )
        {
            glRotated( 180, 1.0, 0., 0. );
        }

        quadObj = gluNewQuadric();
        gluQuadricDrawStyle( quadObj, GLU_FILL );
        gluQuadricNormals( quadObj, GLU_SMOOTH );
        gluCylinder( quadObj, radius, radius, L, 32, 1 );
        gluDeleteQuadric(quadObj);

    glPopMatrix ();
}
