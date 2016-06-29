#ifndef EDGEFACECOLLISION_HH_
#define EDGEFACECOLLISION_HH_

#include "Collision.h"
#include "ElementProxy.h"
class TwistEdgeHandler;
class EdgeFaceCollision: public EdgeCollision, public FaceCollision
{
public:
    EdgeFaceCollision( EdgeProxy* firstProxy, const FaceProxy* face,
            int firstIdx, int secondIdx, int firstApex, int secondApex ) :
            EdgeCollision( firstProxy ), //
            m_mesh( face->getMesh() ), //
            m_face( face ), //
            m_firstIdx( firstIdx ), //
            m_secondIdx( secondIdx ),m_firstApex( firstApex ), m_secondApex( secondApex ),
            m_onBoundary( face->getFace().edgeOnBoundary( m_firstApex ) )
    {

        // Maintain canonical order for sorting
        if ( m_secondIdx < m_firstIdx )
        {
            std::swap( m_secondIdx, m_firstIdx );
            std::swap( m_firstApex, m_secondApex );
        }
    }

    virtual bool analyse();
    virtual bool analyse( TwistEdgeHandler* teh );

    friend bool compare( const EdgeFaceCollision* ef1, const EdgeFaceCollision* ef2 );

    Vec3 meshVelocity( Scalar dt ) const
    {
        return m_meshDisplacement / dt;
    }

    uint32_t faceEdgeId() const
    {
        return m_face->uniqueId() + ( 4 - m_firstApex - m_secondApex ) ;
    }

    const FaceProxy* face() const
    {
        return m_face;
    }

protected:
    void print( std::ostream& os ) const;

    const TriMesh* m_mesh;
    const FaceProxy* m_face;
    int m_firstIdx, m_secondIdx, m_firstApex, m_secondApex;
    bool m_onBoundary ;
    Vec3 m_meshDisplacement;

};

#endif
