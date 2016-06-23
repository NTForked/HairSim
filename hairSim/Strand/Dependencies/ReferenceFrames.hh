/**
 * @file ReferenceFrames.hh
 *
 *  Created on: 17/08/2012
 * \author Jean-Marie Aubry <jaubry@wetafx.co.nz>
 * \note Based on original BASim code from Columbia University
 */

#ifndef REFERENCEFRAMES_HH_
#define REFERENCEFRAMES_HH_

#include "DegreesOfFreedom.hh"

/**
 * \brief This class maintains the first reference frame vectors, orthogonal to the edges.
 *
 * Initialized by space-parallel transportation along the rod, the reference frame vectors later
 * evolve by time-parallel transportation.
 *
 * Unit: no dimension
 */
class ReferenceFrames1: public DependencyNode<Vec3xArray>
{
public:
    ReferenceFrames1( Tangents& tangents ) :
            DependencyNode<Vec3xArray>( 0, tangents.size() ), m_tangents( tangents )
    {
        m_tangents.addDependent( this );

        storeInitialFrames( Vec3x() );
    }

    virtual const char* name() const
    {
        return "ReferenceFrames1";
    }

    /**
     * \brief Propagates the initial (orthornormalized) reference frame initRefFrame1 by
     * space-parallel transportation along the rod
     * @param initRefFrame1 If zero or too far from orthogonal to the first edge, is replaced
     * by an arbitrary orthonormal vector.
     */
    void storeInitialFrames( const Vec3x& initRefFrame1 );

    /**
     * \brief Access to m_previousTangents.
     *
     * m_previousTangents is an internal parameter caching the tangents from the previous
     * time step so we are able to do time-parallel transportation. The only need we have to
     * access it directly is for serialization/deserialization.
     */
    Vec3xArray& getPreviousTangents()
    {
        return m_previousTangents;
    }

    void setPreviousTangents( const Vec3xArray& deserialPrevTangents )
    {
        storeInitialFrames( Vec3x() );
        m_previousTangents = deserialPrevTangents;
    }

    bool checkNormality();

protected:
    /**
     * \brief Computes new reference frames by time-parallel transportation along the
     * m_previousTangents->m_tangents motion.
     */
    virtual void compute();

    Tangents& m_tangents;
    Vec3xArray m_previousTangents;
};

/**
 * \brief The second reference frame is simply maintained as tangent x referenceFrames1
 */
class ReferenceFrames2: public DependencyNode<Vec3xArray>
{
public:
    ReferenceFrames2( Tangents& tangents, ReferenceFrames1& referenceFrames1 ) :
            DependencyNode<Vec3xArray>( 0, tangents.size() ), m_tangents( tangents ), m_referenceFrames1(
                    referenceFrames1 )
    {
        m_tangents.addDependent( this );
        m_referenceFrames1.addDependent( this );
    }

    virtual const char* name() const
    {
        return "ReferenceFrames2";
    }

protected:
    virtual void compute();

    Tangents& m_tangents;
    ReferenceFrames1& m_referenceFrames1;
};

/**
 * \brief This maintains the reference twist defined at the angle between the (space-parallel
 * transported) Bishop frame and the actual (time-parallel transported) reference frame.
 */
class ReferenceTwists: public DependencyNode<std::vector<Scalar> >
{
public:
    ReferenceTwists( Tangents& tangents, ReferenceFrames1& referenceFrames1 ) :
            DependencyNode<std::vector<Scalar> >( 1, tangents.size() ), m_tangents( tangents ), m_referenceFrames1(
                    referenceFrames1 )
    {
        m_tangents.addDependent( this );
        m_referenceFrames1.addDependent( this );
    }

    virtual const char* name() const
    {
        return "ReferenceTwists";
    }

protected:
    virtual void compute();

    Tangents& m_tangents;
    ReferenceFrames1& m_referenceFrames1;
};

#endif /* REFERENCEFRAMES_HH_ */
