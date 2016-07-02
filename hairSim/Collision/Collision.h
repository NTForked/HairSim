#ifndef CONTINUOUS_TIME_COLLISION_H_
#define CONTINUOUS_TIME_COLLISION_H_

#include "../Utils/Definitions.h"

class EdgeProxy;
class FaceProxy;
class TriMesh;
class TwistEdgeHandler;
typedef SparseRowMatx DeformationGradient;

class Collision
{
public:

    Collision( EdgeProxy* edgeProxy );

    virtual ~Collision();

    friend std::ostream& operator<<( std::ostream& os, const Collision& col );

    friend bool compareTimes( const Collision* ct1, const Collision* ct2 );

    EdgeProxy* getFirstEdgeProxy() const
    { return m_firstEdgeProxy; }
    
    Scalar time() const 
    { return m_time; }

    const Vec3& normal() const 
    { return m_normal; }

protected:

    virtual void print( std::ostream& os ) const;

    EdgeProxy* m_firstEdgeProxy;

    Scalar m_time;
    Vec3 m_normal;
};

class EdgeCollision: public Collision
{
public:
    EdgeCollision( EdgeProxy* edgeProxy ):
        Collision( edgeProxy )
    {}

    virtual bool analyse( TwistEdgeHandler* teh ) = 0;

    Scalar abscissa() const
    { return m_s; }

protected:
    Scalar m_s;
};

class FaceCollision
{
public:
    virtual ~FaceCollision()
    {}

    virtual const FaceProxy* face() const = 0;
    Scalar faceFrictionCoefficient() const;

    virtual bool analyse() = 0;

    virtual Vec3 offset( const Vec3& normal ) const;

    Scalar m_u, m_v, m_w;
    Vec3 m_offset;
};

struct CollidingPair
{
public:

    struct Object
    {
        int globalIndex;
        int vertex; // needed for sorting
        Scalar abscissa;
        DeformationGradient* defGrad;
        Vec3 worldVel;
    };

    CollidingPair()
    {}
    
    Scalar m_mu;
    Vec3 m_normal;
    Mat3x m_transformationMatrix;
    std::pair<Object, Object> objects;

    void generateTransformationMatrix();
    void swapIfNecessary();
    bool operator<( const CollidingPair& rhs ) const;
};

typedef std::vector< CollidingPair > CollidingPairs;

bool compareTimes( const Collision* ct1, const Collision* ct2 );


#endif 
