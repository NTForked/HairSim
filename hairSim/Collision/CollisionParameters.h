#ifndef COLLISIONPARAMETERS_HH_
#define COLLISIONPARAMETERS_HH_

#include "../Strand/ElasticStrandParameters.h"

struct CollisionParameters
{
    enum CollisionType
    {
        ROD_ROD, ROD_OTHER
    };

    // For now nothing depends on edges, but maybe it will later

    CollisionParameters() :
        m_useInstantenousForCT( true ),
        m_rejectSelfCollisions( true ), 
        m_collisionRadius( .1 ), 
        m_externalCollisionsRadius( .1 ), 
        m_frictionCoefficient( .2 ), 
        m_rootImmunityLength( 0.0 ), 
        m_associatedStrandParams( NULL )
    {}

    CollisionParameters( 
        Scalar selfCollisionRadius, 
        Scalar externalCollisionsRadius,
        Scalar frictionCoefficient, 
        Scalar rootCollisionImmunity, 
        bool rejectSelfCollisions = true ):
            m_useInstantenousForCT( true ),
            m_rejectSelfCollisions( rejectSelfCollisions ), 
            m_collisionRadius( selfCollisionRadius ), 
            m_externalCollisionsRadius( externalCollisionsRadius ), 
            m_frictionCoefficient( frictionCoefficient ),
            m_rootImmunityLength( rootCollisionImmunity ), 
            m_associatedStrandParams( NULL )
    {}

    Scalar externalCollisionsRadius( unsigned edgeIdx, Scalar angle = 0. ) const
    {
        return collisionRadius( ROD_OTHER, edgeIdx, angle );
    }

    Scalar collisionRadius( unsigned edgeIdx, Scalar angle = 0. ) const
    {
        return collisionRadius( ROD_ROD, edgeIdx, angle );
    }

    Scalar collisionRadius( CollisionType type, unsigned edgeIdx, Scalar angle = 0. ) const
    {
        const Scalar baseRadius = type == ROD_ROD ? m_collisionRadius : m_externalCollisionsRadius;

        if ( m_associatedStrandParams )
        {
            const Scalar majorRad = m_associatedStrandParams->interpolatedRadiusMultiplier( edgeIdx ) * baseRadius;
            return majorRad;
        }
        return baseRadius;
    }

    Scalar getLargerCollisionsRadius( unsigned edgeIdx, Scalar angle = 0. ) const
    {
        Scalar rod = collisionRadius( ROD_ROD, edgeIdx, angle );
        Scalar other = collisionRadius( ROD_OTHER, edgeIdx, angle );
        return rod > other ? rod : other;
    }

    void multiplyExternalCollisionsRadius( Scalar m )
    {
        this->m_externalCollisionsRadius *= m;
    }

    void multiplyFrictionCoefficient( Scalar m )
    {
        this->m_frictionCoefficient *= m;
    }

    void multiplyRootImmunityLength( Scalar m )
    {
        this->m_rootImmunityLength *= m;
    }

    void multiplyCollisionRadius( Scalar m )
    {
        this->m_collisionRadius *= m;
    }

    void setAssociatedStrandParameters( const ElasticStrandParameters& strandParams )
    {
        m_associatedStrandParams = &strandParams;
    }

//////////

    bool m_useInstantenousForCT; // still use CT's BVH culling, but Instantenous instead of CT in narrow phase
    bool m_rejectSelfCollisions;

// private:

    Scalar m_collisionRadius; // this rod's radius against other rods
    Scalar m_externalCollisionsRadius; // this rod's radius against other collisions (mesh)
    Scalar m_frictionCoefficient;
    Scalar m_rootImmunityLength;

    // Required for getting radius interpolation
    const ElasticStrandParameters* m_associatedStrandParams;
};

#endif /* COLLISIONPARAMETERS_HH_ */
