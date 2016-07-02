#ifndef VERTEXFACECOLLISION_HH_
#define VERTEXFACECOLLISION_HH_

#include "Collision.h"
#include "ElementProxy.h"

class VertexFaceCollision: public Collision, public FaceCollision
{
public:
    VertexFaceCollision( EdgeProxy* firstProxy, const FaceProxy* faceProxy ) :
            Collision( firstProxy ), 
            m_faceProxy( faceProxy )
    {
        m_firstStrand = firstProxy->getStrandPointer();
        m_firstVertex = firstProxy->getVertexIndex();

    }

    virtual ~VertexFaceCollision()
    {}

    virtual bool analyse();

    friend bool compare( const VertexFaceCollision* vf1, const VertexFaceCollision* vf2 );

    Vec3 meshVelocity( Scalar dt ) const
    { return m_meshDisplacement / dt; }

    const FaceProxy* face() const
    { return m_faceProxy; }

protected:
    void print( std::ostream& os ) const;

    const FaceProxy* const m_faceProxy;
    Vec3 m_meshDisplacement;
    Vec3 m_collisionOffset;

    ElasticStrand* m_firstStrand;
    int m_firstVertex;
};

#endif
