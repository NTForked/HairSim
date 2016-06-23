#ifndef CAMERA_HH
#define CAMERA_HH

#include "../../Utils/TypeDefs.h"

/** Class that holds camera parameters written by Adrian Secord, with
    modifications by Miklos Bergou. */
class Camera
{
public:

  /// Projection mode options
  enum ProjMode {
    ORTHOGRAPHIC,                               ///< Orthogonal projection
    PERSPECTIVE                                 ///< Perspective projection
  };

  /// Default constructor
  Camera();

  /// Copy constructor
  Camera(const Camera& other);

  /// Assignment operator
  Camera& operator= (const Camera& other);

  /// Destructor
  ~Camera() {}

  /// Set the viewport to an orthographic view.
  void setOrthographic(const scalar left, const scalar right,
                       const scalar top, const scalar bottom);

  /// Set the viewport to a perspective view.
  void setPerspective(const scalar fovy, const scalar aspect);

  /// Get the view center position
  const Vec3& getViewCenter() const;

  /// Get the eye position
  const Vec3& getEye() const;

  /// Get the up vector
  const Vec3& getUp() const;

  /// Set the view center position
  void setViewCenter(const Vec3& v);

  /// Set the eye position
  void setEye(const Vec3& e);

  /// Set the up vector
  void setUp(const Vec3& u);

  /// Get the view parameters.
  void getViewParams(ProjMode* m, std::vector< scalar >* p) const;

  /// Set the viewport
  void setViewport(const int width, const int height);

  /// Get the viewport
  void getViewport(int* width, int* height) const;

  /// Set the z-clipping planes.
  void setZClipping( const scalar near, const scalar far);

  /// Get the z-clipping planes.
  void getZClipping( scalar& near, scalar& far) const;

  /// Apply the current camera projection to the OpenGL state.
  void applyProjection() const;

  /// Apply the viewport to the OpenGL state.
  void applyViewport() const;

  /// Apply the current camera viewing position/direction to the OpenGL state.
  void applyCamera() const;

  /** \name Camera positioning and movements. */
  //@{

  /// Set a default view for a 3D scene.
  /// \c sceneRadius should be the bounding radius of the scene.
  void setDefault3D(const scalar sceneRadius);

  /// Set a default view for a 2D scene in the z == 0 plane.
  /// \c sceneRadius should be the bounding radius of the scene.
  void setDefault2D(const scalar sceneRadius);

  /// Translate the camera by a vector.
  void translateEye(const Vec3& v);

  /// Translate the view center by a vector.
  void translateCenter(const Vec3& v);

  /// Rotate the camera about the view center by a 4x4 rotation matrix.
  void orbit(const scalar m[4][4]);

  //@}

  /// Get a set of vectors that span the view plane in world coords.
  /// The first vector points "right", the second points "up", and the
  /// third vector points towards the view center.
  void getSpanningSet( Vec3* right, Vec3* up, Vec3* center) const;

protected:

  Vec3 m_viewCenter;                    ///< World-space point to rotate about
  Vec3 m_eye;                           ///< World-space location of camera
  Vec3 m_up;                            ///< Direction of the top of camera
  Vec3 m_dir;                           ///< Direction the camera is looking at

  ProjMode m_projMode;                   ///< Projection mode
  std::vector< scalar > m_projParams;      ///< Projection parameters
  std::vector< scalar> m_zClipping;       ///< Near and far clipping planes
  std::vector< int > m_viewport;           ///< Viewport size
};

#endif // CAMERA_HH
