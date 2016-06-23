
#ifndef TWIST_EDGE_COLLISION_HH_
#define TWIST_EDGE_COLLISION_HH_

#include "EdgeEdgeCollision.hh"
#include "ElementProxy.hh"
#include "CTCD.h"
#include "../Dynamic/ImplicitStepper.hh"
#include "TwistEdgeHandler.hh"

namespace strandsim
{

class TwistEdgeCollision: public EdgeEdgeCollision
{
public:
    TwistEdgeCollision( ElasticStrand* firstStrand, int firstVertex, ElasticStrand* secondStrand,
            int secondVertex, ImplicitStepper* firstStepper, ImplicitStepper* secondStepper,
            TwistEdge* firstProxy, TwistEdge* secondProxy )
    : 
        EdgeEdgeCollision( firstStrand, firstVertex, secondStrand, secondVertex, firstStepper, secondStepper ),
        m_firstProxy( firstProxy ),
        m_secondProxy( secondProxy )
    {
        // Maintain canonical order for sorting
        if ( m_firstStrand != firstStrand ){
            std::swap( m_firstProxy, m_secondProxy );
        }
    }

    bool partnerParentsColliding( TwistEdgeHandler* teh );

    virtual bool analyse( TwistEdgeHandler* teh );

    bool operator <( EdgeEdgeCollision* other) const{
        return compare( this, other);
    }
    
    TwistEdge* m_firstProxy;
    TwistEdge* m_secondProxy;
    Scalar twistAngle;
};

}

#endif
