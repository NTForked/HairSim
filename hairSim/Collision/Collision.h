#ifndef CONTINUOUS_TIME_COLLISION_H_
#define CONTINUOUS_TIME_COLLISION_H_

class EdgeProxy;
class FaceProxy;
class TriangularMesh;
typedef SparseRowMatx DeformationGradient;

bool compareTimes( const Collision* ct1, const Collision* ct2 );

class Collision
{
public:

    Collision( EdgeProxy* edgeProxy );

    virtual ~Collision();

    friend std::ostream& operator<<( std::ostream& os, const Collision& col );

    virtual bool analyse() = 0;

    friend bool compareTimes( const CollisionBase* ct1, const CollisionBase* ct2 );

    EdgeProxy* getFirstEdgeProxy() const
    { return m_firstEdgeProxy; }
    
    Scalar time() const 
    { return m_time; }

    const Vec3x& normal() const 
    { return m_normal; }

protected:

    virtual void print( std::ostream& os ) const;

    EdgeProxy* m_firstEdgeProxy;

    Scalar m_time;
    Vec3x m_normal;
};

class EdgeCollision: public Collision
{
public:
    EdgeCollision( EdgeProxy* edgeProxy ):
        Collision( edgeProxy )
    {}

    Scalar abscissa() const
    { return m_s; }

private:
    Scalar m_s;
};

class FaceCollision
{
public:
    virtual ~FaceCollision()
    {}

    virtual const FaceProxy* face() const = 0;
    Scalar faceFrictionCoefficient() const;

    virtual Vec3x offset() const;

    Scalar m_u, m_v, m_w;
    Vec3x m_offset;
};

struct CollisionPair
{
public:

    struct Object
    {
        DeformationGradient* defGrad;
        Vec3x worldVel;
    };

    ProximityCollision():
    {}
    
    Vec3x m_warmStartImpulse;
    Mat3x m_transformationMatrix;

    std::pair<Object, Object> objects;

    void generateTransformationMatrix();
  
    void print( std::ostream& os ) const;
    void printShort( std::ostream& os ) const;
    bool operator< ( const ProximityCollision &rhs ) const;
    void swapIfNecessary();
};

typedef std::vector<ProximityCollision> ProximityCollisions;


#endif 
