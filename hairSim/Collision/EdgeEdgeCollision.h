#ifndef EDGE_EDGE_COLLISION_H_
#define EDGE_EDGE_COLLISION_H_

#include "Collision.h"

class TwistEdgeHandler;

class EdgeEdgeCollision: public EdgeCollision
{
public:
    EdgeEdgeCollision( EdgeProxy* firstEdgeP, EdgeProxy* secondEdgeP ):
            EdgeCollision( firstEdgeP ), 
            m_secondEdgeProxy( secondEdgeP )
    {
        if ( m_secondEdgeProxy < m_firstEdgeProxy )
        { // Maintain canonical order for sorting
            std::swap( m_firstEdgeProxy, m_secondEdgeProxy );
        }
    }

    virtual bool analyse( TwistEdgeHandler* teh );

    bool operator <( EdgeEdgeCollision* other ) const
    { return compare( this, other); }

    friend bool compare( const EdgeEdgeCollision* ee1, const EdgeEdgeCollision* ee2 );
    
    EdgeProxy* getSecondEdgeProxy() const
    { return m_secondEdgeProxy; }

    Scalar getFirstAbscissa() const
    { return m_s; }
    
    Scalar getSecondAbscissa() const
    { return m_t; }

protected:

    void validCollision( TwistEdgeHandler* teh );
    void print( std::ostream& os ) const;
    void printShort( std::ostream& os ) const;

    EdgeProxy* m_secondEdgeProxy;
    Scalar m_t;
};

#endif
