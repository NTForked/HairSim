#ifdef __APPLE__
#  include <OpenGL/gl.h>
#  include <OpenGL/glu.h>
#  include <GLUT/glut.h>
#else
#  include <GL/gl.h>
#  include <GL/glu.h>
#  include <GL/glut.h>
#endif

#include "ViewController.h"

const Scalar ViewController::eps = std::numeric_limits<Scalar>::epsilon();

ViewController::ViewController()
  : m_centerMode(CENTER_OBJECT)
  , m_trackball(&m_camera)
  , m_translator(&m_camera)
  , m_zoomer(&m_camera)
{
  m_camera.setPerspective(60, 1);
  m_camera.setViewport(100, 100);
}

void ViewController::ApplyCamera()
{
  m_camera.applyViewport();
  m_camera.applyProjection();

  glMatrixMode( GL_MODELVIEW );
  glLoadIdentity();

  m_camera.applyCamera();
}

void ViewController::setCamera(const Camera& c) {
  m_camera = c;
}

const Camera& ViewController::getCamera() const {
  return m_camera;
}

Camera& ViewController::getCamera() {
  return m_camera;
}

void ViewController::setCenterMode(const CenterMode m) {
  m_centerMode = m;
  switch (m_centerMode) {
  case CENTER_WORLD_ORIGIN:
    m_camera.setViewCenter(Vec3(0,0,0));
    break;

  case CENTER_OBJECT:
    m_camera.setViewCenter(m_objCenter);
    break;
  }
}

void ViewController::setViewCenter(const Vec3& p) {
  m_objCenter = p;
  m_camera.setViewCenter(p);
}

void ViewController::getCenterOfInterest(ViewController::CenterMode* m,
                                         Vec3* p) const
{
  if (m)
    *m = m_centerMode;

  if (p)
    *p = m_camera.getViewCenter();
}

void ViewController::setViewDirection(const Vec3& d) {
  m_camera.setEye(m_camera.getViewCenter() - d);
  m_camera.setUp(Vec3(0, 1, 0));
}

void ViewController::setBoundingRadius(Scalar r) {
  m_boundingRadius = r;
  m_camera.setDefault3D(m_boundingRadius);
  m_translator.setScale(2 * m_boundingRadius);
  m_zoomer.setScale(2 * m_boundingRadius);
}

Scalar ViewController::getBoundingRadius() {
  return m_boundingRadius;
}

void ViewController::beginRotationDrag(const Scalar x, const Scalar y) {
  m_trackball.start(Vec2(x,y));
}

void ViewController::endRotationDrag(const Scalar x, const Scalar y) {
  m_trackball.update(Vec2(x,y));
  m_trackball.stop();
}

void ViewController::beginTranslationDrag(const Scalar x, const Scalar y) {
  m_translator.start(Vec2(x,y));
}

void ViewController::endTranslationDrag(const Scalar x, const Scalar y) {
  m_translator.update(Vec2(x,y));
  m_translator.stop();
}

void ViewController::beginZoomDrag(const Scalar x, const Scalar y) {
  m_zoomer.start(Vec2(x,y));
}

void ViewController::endZoomDrag(const Scalar x, const Scalar y) {
  m_zoomer.update(Vec2(x,y));
  m_zoomer.stop();
}

// beginScriptingDrag
// endScriptingDrag
// add to updateDrag

void ViewController::updateDrag(const Scalar x, const Scalar y) {
  const Vec2 v(x,y);
  m_trackball.update(v);
  m_translator.update(v);
  m_zoomer.update(v);
}

void ViewController::calcTranslation(const Vec2& start, const Vec2& stop,
                                     const Vec3& viewCenter, const Vec3& eye,
                                     const Vec3& up, const Scalar scale,
                                     Vec3& translation)
{
  // Points to the object
  const Vec3 n = (viewCenter - eye).normalized();
  // Spans the view plane in world coords
  const Vec3 w2 = (up - up.dot(n) * n).normalized();
  const Vec3 w1 = n.cross(w2);

  const Vec2 v = (stop - start).normalized() * scale;

  translation = v[0] * w1 + v[1] * w2;

  assert(translation.dot(n) < 10*eps);
}

