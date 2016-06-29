#include "Translator.h"

Translator::Translator(Camera* c, const Scalar s)
  : m_camera(c)
  , m_translating(false)
  , m_translation(0,0,0)
  , m_scale(s)
{}

void Translator::setCamera(Camera* c)
{
  m_camera = c;
}

void Translator::setScale(const Scalar s)
{
  m_scale = s;
}

void Translator::start(const Vec2& p)
{
  m_translating = true;
  m_startPos = p;
}

void Translator::update(const Vec2& p)
{
  if (!m_translating)
    return;

  Vec3 right, up;
  m_camera->getSpanningSet(&right, &up, NULL);

  const Vec2 v = (p - m_startPos) * m_scale;
  m_translation = v[0] * right + v[1] * up;

  m_camera->translateEye(-m_translation);
  m_camera->translateCenter(-m_translation);

  m_startPos = p;
}

void Translator::stop()
{
  m_translating = false;
}

const Vec3& Translator::getTranslation() const
{
  return m_translation;
}

