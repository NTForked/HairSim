/*
 * GravitationForce.hh
 *
 *  Created on: 14/07/2011
 *      Author: Jean-Marie Aubry <jaubry@wetafx.co.nz>
 */

#ifndef GRAVITATIONFORCE_HH_
#define GRAVITATIONFORCE_HH_

#include "ForceBase.hh"
#include "../Utils/Definitions.h"
#include "../Strand/ElasticStrand.h"
 
class GravitationForce: public ForceBase
{
public:
    static const IndexType s_first = 0; // The first index on which this force can apply
    static const IndexType s_last = 0; // The last index (counting from the end)

    typedef Vec3 LocalForceType;
    typedef Mat3x LocalJacobianType;

    GravitationForce();
    virtual ~GravitationForce();

    static std::string getName()
    {
        return "gravitation";
    }

    static Scalar localEnergy( const ElasticStrand& strand, const StrandState& geometry,
            const IndexType vtx );

    template<typename LocalT>
    static void computeLocal( LocalT& local, const ElasticStrand& strand,
            const StrandState& geometry, const IndexType vtx );

    template<typename GlobalT, typename LocalT>
    static void addInPosition( GlobalT& global, const IndexType vtx, const LocalT& local );

    static void setGravity( const Vec3& gravity )
    {
        s_gravity = gravity;
    }

private:
    static Vec3 s_gravity;
};

#endif /* GRAVITATIONFORCE_HH_ */
