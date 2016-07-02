#ifndef RENDER_UTILS_H
#define RENDER_UTILS_H

#include "../Utils/Definitions.h"
#include "OpenGLDecl.h"

void renderSphere( const Vec3& p, GLdouble radius );

void renderCylinder( const Vec3& p, const Vec3& q, GLdouble radius );

#endif
