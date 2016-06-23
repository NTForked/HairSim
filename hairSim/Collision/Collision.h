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

    virtual Vec3x offset() const;

    EdgeProxy* getFirstEdgeProxy() const
    { return m_firstEdgeProxy; }
    
    Scalar time() const 
    { return m_time; }

    const Vec3x &normal() const 
    { return m_normal; }

protected:

    virtual void print( std::ostream& os ) const;

    EdgeProxy* m_firstEdgeProxy;

    Scalar m_time;
    Scalar m_mu;
    Vec3x m_normal;
    Vec3x m_offset;
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

    Scalar m_u, m_v, m_w;
};


struct CollisionPair
{
public:

    struct Object
    {

        DeformationGradient* defGrad;
        Vec3x freeVel;
    };

    ProximityCollision():
    {}
    
    Vec3x warmStartImpulse;
    Mat3x m_transformationMatrix;

    std::pair<Object, Object> objects;

    void generateTransformationMatrix() ;
    void updateTransformationMatrix( const Mat3x& previous ) ;
  
    void print( std::ostream& os ) const;
    void printShort( std::ostream& os ) const;
    bool operator< ( const ProximityCollision &rhs ) const;
    void swapIfNecessary();
};

typedef std::vector<ProximityCollision> ProximityCollisions;


#endif 
