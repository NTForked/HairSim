#include "TrackBall.h"

namespace {

inline void
rotate(Scalar m[4][4], const Vec3& axis,
       const Scalar radians)
{
  Vec3 naxis( axis );
  naxis.normalize();

  Scalar x = naxis[0];
  Scalar y = naxis[1];
  Scalar z = naxis[2];

  Scalar c = cos( radians );
  Scalar s = sin( radians );
  Scalar t = 1 - c;

  m[0][0] = t*x*x + c;
  m[0][1] = t*x*y - z*s;
  m[0][2] = t*x*z + y*s;
  m[1][0] = t*x*y + z*s;
  m[1][1] = t*y*y + c;
  m[1][2] = t*y*z - x*s;
  m[2][0] = t*x*z - y*s;
  m[2][1] = t*y*z + x*s;
  m[2][2] = t*z*z + c;

  m[0][3] = 0;
  m[1][3] = 0;
  m[2][3] = 0;
  m[3][0] = 0;
  m[3][1] = 0;
  m[3][2] = 0;
  m[3][3] = 1;
}

} // namespace (global)

TrackBall::TrackBall(Camera* c)
  : m_camera(c)
  , m_rotating(false)
  , m_mode(GIMBAL)
{}

void TrackBall::start(const Vec2& p)
{
  assert(p[0] >= -1 && p[0] <= 1 && p[1] >= -1 && p[1] <= 1);
  m_rotating = true;
  m_startPos = p;
  for (int i = 0; i < 4; ++i)
    m_rotation[i] = 0;
}

void TrackBall::update(const Vec2& p)
{
  if (!m_rotating)
    return;

  Scalar m[4][4];
  if (m_startPos != p) {
    const Scalar coef(M_PI/2.0);

    const Scalar left_right_motion = p[0] - m_startPos[0];
    const Scalar up_down_motion = p[1] - m_startPos[1];

    // rotate about the 'up' vector
    Vec3 up;
    if (m_mode == GIMBAL) up = Vec3(0,1,0);
    else                  up = m_camera->getUp();
    int sign = (up.dot(m_camera->getUp()) > 0) ? 1 : -1;
    rotate(m, up, coef * sign * left_right_motion);
    m_camera->orbit(m);

    // rotate about the horizontal vector
    Vec3 horizontal =
      m_camera->getUp().cross(m_camera->getEye() - m_camera->getViewCenter());
    rotate(m, horizontal, -coef * up_down_motion);
    m_camera->orbit(m);

    m_startPos = p;
  }
}

void TrackBall::stop()
{
  m_rotating = false;
}

void TrackBall::getRotation(Scalar r[4]) const
{
  for (int i = 0; i < 4; ++i) {
    r[i] = m_rotation[i];
  }
}

TrackBall::TrackBallMode TrackBall::getTrackBallMode() const
{
  return m_mode;
}

void TrackBall::setTrackBallMode(const TrackBall::TrackBallMode m)
{
  m_mode = m;
}
