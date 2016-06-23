#ifndef VIEWCONTROLLER_HH
#define VIEWCONTROLLER_HH

#include "Camera.h"
#include "TrackBall.h"
#include "Translator.h"
#include "Zoomer.h"

class ViewController
{
public:

  /// Center of interest options
  enum CenterMode {
    CENTER_WORLD_ORIGIN,                   ///< Center the view on the origin
    CENTER_OBJECT,                         ///< Center the view on the mesh
  };

  /// Default constructor.
  ViewController();

  void ApplyCamera();

  /// Set the current camera.
  void setCamera(const Camera& c);

  /// Get the current camera.
  const Camera& getCamera() const;

  /// Get the current camera.
  Camera& getCamera();

  /// Set the view center of interest
  void setCenterMode(const CenterMode m);
  void setViewCenter(const Vec3& p);

  void setBoundingRadius(scalar r);
  scalar getBoundingRadius();

  /// Get the view center of interest.  Either parameter can be NULL if
  /// it is not needed.
  void getCenterOfInterest(CenterMode* m, Vec3* p = NULL) const;

  /// Set a particular viewing direction
  void setViewDirection(const Vec3& d);

  /** \name Interaction commands

      Commands to modify the camera from user interaction. */

  //@{

  /// Begin a rotation drag at a point in pixel coordinates.
  void beginRotationDrag(const scalar x, const scalar y);

  /// End a rotation drag at a point in pixel coordinates.
  void endRotationDrag(const scalar x, const scalar y);

  /// Begin a translation drag at a point in pixel coordinates.
  void beginTranslationDrag(const scalar x, const scalar y);

  /// End a translation drag at a point in pixel coordinates.
  void endTranslationDrag(const scalar x, const scalar y);

  /// Begin a zoom drag at a point in pixel coordinates.
  void beginZoomDrag(const scalar x, const scalar y);

  /// End a zoom drag at a point in pixel coordinates.
  void endZoomDrag(const scalar x, const scalar y);

  /// Update a current drag operation with new coordinates.
  void updateDrag(const scalar x, const scalar y);

  //@}

private:

  /// Smallest amount that can be added to a scalar and cause a difference.
  static const scalar eps;

  /// Compute the translation needed to move a point from world-space start to
  /// intersect the perpendicular line that pierces the screen at (x,y).
  static void calcTranslation(const Vec2& start,
                              const Vec2& stop,
                              const Vec3& viewCenter,
                              const Vec3& eye,
                              const Vec3& up,
                              const scalar scale,
                              Vec3& translation);

  /// Not copyable
  ViewController(const ViewController& m);

  /// Not copyable
  ViewController& operator= (const ViewController& m);

  Camera m_camera;                       ///< Camera for viewing scene.

  CenterMode m_centerMode;               ///< Centering mode

  Vec3 m_objCenter;                     ///< Center of the object.
  scalar m_boundingRadius;               ///< Bounding radius of mesh.

  TrackBall m_trackball;                 ///< Trackball for user rotation.
  Translator m_translator;               ///< User translations.
  Zoomer m_zoomer;                       ///< User zooms.

};

#endif // VIEWCONTROLLER_HH
