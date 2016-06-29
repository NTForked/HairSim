#include "Zoomer.h"

Zoomer::Zoomer(Camera* c, const Scalar s)
  : m_camera(c)
  , m_translating(false)
  , m_scale(s)
{}

void Zoomer::setCamera(Camera* c)
{
  m_camera = c;
}

void Zoomer::setScale(const Scalar s)
{
  m_scale = s;
}

void Zoomer::start(const Vec2& p)
{
  m_translating = true;
  m_startPos = p;
}

void Zoomer::update(const Vec2& p)
{
  if (!m_translating) return;

  assert(m_camera);
  Vec3 in;
  m_camera->getSpanningSet(NULL, NULL, &in);

  const Vec3 translation = in * m_scale * (p[1] - m_startPos[1]);

  m_camera->translateEye(translation);

  m_startPos = p;
}

void Zoomer::stop()
{
  m_translating = false;
}
