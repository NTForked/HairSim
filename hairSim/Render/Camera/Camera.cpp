#ifdef __APPLE__
#  include <OpenGL/gl.h>
#  include <OpenGL/glu.h>
#  include <GLUT/glut.h>
#else
#  include <GL/gl.h>
#  include <GL/glu.h>
#  include <GL/glut.h>
#endif

#include "Camera.h"

Camera::Camera()
  : m_viewCenter(0,0,-1)
  , m_eye(0,0,0)
  , m_up(0,1,0)
  , m_dir(0,0,-1)
  , m_projMode(Camera::PERSPECTIVE)
  , m_projParams(2)
  , m_zClipping(2)
{
  m_projParams[0] = 54.43;
  m_projParams[1] = 1;
  m_zClipping[0] = 0.1;
  m_zClipping[1] = 1000;
}

Camera::Camera(const Camera& other)
  : m_viewCenter(other.m_viewCenter)
  , m_eye(other.m_eye)
  , m_up(other.m_up)
  , m_dir(other.m_dir)
  , m_projMode(other.m_projMode)
  , m_projParams(other.m_projParams)
  , m_zClipping(other.m_zClipping)
{}

Camera& Camera::operator=(const Camera& other)
{
  m_viewCenter = other.m_viewCenter;
  m_eye = other.m_eye;
  m_up = other.m_up;
  m_dir = other.m_dir;
  m_projMode = other.m_projMode;
  m_projParams = other.m_projParams;
  m_zClipping = other.m_zClipping;
  return *this;
}

void Camera::setOrthographic(const scalar left, const scalar right,
                             const scalar top, const scalar bottom)
{
  m_projMode = ORTHOGRAPHIC;
  m_projParams.resize(4);
  m_projParams[0] = left;
  m_projParams[1] = right;
  m_projParams[2] = top;
  m_projParams[3] = bottom;
}

void Camera::setPerspective(const scalar fovy, const scalar aspect)
{
  m_projMode = PERSPECTIVE;
  m_projParams.resize(2);
  m_projParams[0] = fovy;
  m_projParams[1] = aspect;
}

const Vec3& Camera::getViewCenter() const
{
  return m_viewCenter;
}

const Vec3& Camera::getEye() const
{
  return m_eye;
}

const Vec3& Camera::getUp() const
{
  return m_up;
}

void Camera::setViewCenter(const Vec3& c)
{
  m_viewCenter = c;
  m_dir = (m_viewCenter - m_eye).normalized();
}

void Camera::setEye(const Vec3& e)
{
  m_eye = e;
}

void Camera::setUp(const Vec3& u)
{
  m_up = u;
}

void Camera::getViewParams(Camera::ProjMode* m, std::vector<scalar>* p) const
{
  assert(m && p);
  *m = m_projMode;
  p->resize(m_projParams.size());
  std::copy(m_projParams.begin(), m_projParams.end(), p->begin());
}

void Camera::setViewport(const int width, const int height)
{
  m_viewport.resize(2);
  m_viewport[0] = width;
  m_viewport[1] = height;
}

void Camera::getViewport(int* width, int* height) const
{
  assert(m_viewport.size() == 2);
  if (width) *width = m_viewport[0];
  if (height) *height = m_viewport[1];
}

void Camera::setZClipping(const scalar near, const scalar far)
{
  assert(m_zClipping.size() == 2);
  m_zClipping[0] = std::min(m_zClipping[0], near);
  m_zClipping[1] = std::max(m_zClipping[1], far);
}

void Camera::getZClipping(scalar& near, scalar& far) const
{
  assert(m_zClipping.size() == 2);
  near = m_zClipping[0];
  far = m_zClipping[1];
}

void Camera::setDefault3D(const scalar sceneRadius)
{
  m_dir = Vec3(-28, -21, -28).normalized();
  m_eye = m_viewCenter - sceneRadius * 3.0f * m_dir;
  //m_eye = m_viewCenter + sceneRadius * 3.0f * Vec3(1,1,1).normalized();
  //m_eye = m_viewCenter + sceneRadius * 3.0f * Vec3(0,0,1).normalized();
  m_up = Vec3(0,1,0);
}

void Camera::setDefault2D(const scalar sceneRadius)
{
  m_eye = m_viewCenter + sceneRadius * 3.0f * Vec3(0,0,1).normalized();
  m_up = Vec3(0,1,0);
}

void Camera::translateEye(const Vec3& v)
{
  m_eye += v;
}

void Camera::translateCenter(const Vec3& v)
{
  m_viewCenter += v;
}

// The lame-o scalar pointer representation matches some legacy code we use.
void Camera::orbit(const scalar m[4][4])
{
  {
    const Vec3 e = m_eye - m_viewCenter;
    m_eye[0] = e[0] * m[0][0] + e[1] * m[1][0] + e[2] * m[2][0];
    m_eye[1] = e[0] * m[0][1] + e[1] * m[1][1] + e[2] * m[2][1];
    m_eye[2] = e[0] * m[0][2] + e[1] * m[1][2] + e[2] * m[2][2];
    m_eye += m_viewCenter;
  }

  {
    const Vec3 u = m_up;
    m_up[0] = u[0] * m[0][0] + u[1] * m[1][0] + u[2] * m[2][0];
    m_up[1] = u[0] * m[0][1] + u[1] * m[1][1] + u[2] * m[2][1];
    m_up[2] = u[0] * m[0][2] + u[1] * m[1][2] + u[2] * m[2][2];
  }

  {
    const Vec3 d = m_dir;
    m_dir[0] = d[0] * m[0][0] + d[1] * m[1][0] + d[2] * m[2][0];
    m_dir[1] = d[0] * m[0][1] + d[1] * m[1][1] + d[2] * m[2][1];
    m_dir[2] = d[0] * m[0][2] + d[1] * m[1][2] + d[2] * m[2][2];
  }
}

void Camera::getSpanningSet(Vec3* r, Vec3* u, Vec3* c) const
{
  //const Vec3 n = (m_viewCenter - m_eye).normalized();
  const Vec3 n = m_dir;
  const Vec3 w2 = (m_up - m_up.dot(n) * n).normalized();
  const Vec3 w1 = n.cross(w2);

  if (r) *r = w1;
  if (u) *u = w2;
  if (c) *c = n;
}

void Camera::applyProjection() const
{
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  if (m_projMode == ORTHOGRAPHIC) {
    assert(m_projParams.size() == 4);
    glOrtho(m_projParams[0], m_projParams[1],
            m_projParams[2], m_projParams[3],
            m_zClipping[0], m_zClipping[1]);

  } else if (m_projMode == PERSPECTIVE) {
    assert(m_projParams.size() == 2);
    gluPerspective(m_projParams[0],
                   m_projParams[1] * m_viewport[0] / (GLfloat)m_viewport[1],
                   m_zClipping[0], m_zClipping[1]);

  } else {
      std::cout << "Camera Error : Unknown mode\n";
  }

  glMatrixMode(GL_MODELVIEW);
}

void Camera::applyViewport() const
{
  glViewport(0, 0, m_viewport[0], m_viewport[1]);
}

void Camera::applyCamera() const
{
  Vec3 center = m_eye + m_dir;
  gluLookAt(m_eye[0], m_eye[1], m_eye[2],                       // Eye position
            center[0], center[1], center[2],              // Look-at position
            m_up[0], m_up[1], m_up[2]);                         // Up vector
}

