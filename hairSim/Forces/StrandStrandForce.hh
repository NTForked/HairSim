/*
 * StrandStrandForce.hh
 *
 *  Created on: 25/10/2011
 *      Author: Jean-Marie Aubry <jaubry@wetafx.co.nz>
 */

#ifndef STRANDSTRANDFORCE_HH_
#define STRANDSTRANDFORCE_HH_

#include "ForceBase.hh"
#include "../Static/StaticCollision.hh"
#include <tr1/memory>

namespace strandsim
{

class StrandStrandForce: public RuntimeCollisionForceBase
{
public:
    // Constructor computing the force orientation from current position
    StrandStrandForce( const ElasticStrand* strand1, int edge1, const ElasticStrand* strand2,
            int edge2, double thickness, double stiffness, double hardness, double friction );

    // Constructor knowing the force orientation
    StrandStrandForce( const ElasticStrand* strand1, int edge1, const ElasticStrand* strand2,
            int edge2, double thickness, double stiffness, double hardness, double friction,
            double initialSign );

    virtual ~StrandStrandForce();

    virtual std::string getName() const;

    virtual void accumulateEFJ( StrandState& geometry, const ElasticStrand& strand );

    virtual void accumulateEF( StrandState& geometry, const ElasticStrand& strand );


    virtual void accumulateJ( StrandState& geometry, const ElasticStrand& strand );

    // New style interface (consistent with static-class forces)
    virtual void accumulateFutureE( Scalar& energy, const ElasticStrand& strand );
    virtual void accumulateCurrentF( VecXx& globalF, const ElasticStrand& strand );


    // RuntimeCollisionForceBase
    virtual Vec3x getWorldCollisionPoint( const StaticCollision &collision,
                                          const ElasticStrand& strand ) const ;
    virtual Vec3x getWorldCollisionNormal( const StaticCollision &collision,
                                           const ElasticStrand& strand ) const ;

    virtual bool collidesWith( const ElasticStrand& strand ) const
    {
        return &strand == m_strandP || &strand == m_strandQ ;
    }

    void updateCollisionData() ;
    void registerCollision( const ElasticStrand& strand, StrandState& geometry ) const;

    bool isActive() const ;

    friend class StrandStrandForceExactComparator;
    friend class StrandStrandForceProximityComparator;

public:
    Scalar tetraVolume() const;
    Vec3x gradP0TetraVolume() const;
    Vec3x gradP1TetraVolume() const;
    Vec3x gradQ0TetraVolume() const;
    Vec3x gradQ1TetraVolume() const;
    Mat3x hessP0P1TetraVolume() const;
    Mat3x hessQ0Q1TetraVolume() const;
    bool SegmentsProjOnEachOther() const;
    std::pair<int, int> SegmentsSituation() const;
    Scalar getInitialSign() const;
    void findTheVertices( const StrandState& geometry ) const;

    int& getEdgeP();
    int& getEdgeQ();
    const ElasticStrand *getStrandP() const;
    const ElasticStrand *getStrandQ() const;
    const ElasticStrand *getOther( const ElasticStrand * me ) const;

private:

    void GetClosestPoints( Scalar &s, Scalar &t ) const;

    typedef enum
    {
        P, Q
    } ClientNumber; // Who the force is actually applied to

    void accumulateEF( ClientNumber clientStrand, VecXx& globalF, Scalar& globalE ) const ;
    void accumulateE( ClientNumber clientStrand,  Scalar& globalE ) const ;

    const ElasticStrand* m_strandP;
    const ElasticStrand* m_strandQ;
    int m_edgeP;
    int m_edgeQ;
    Scalar m_initialSign;
    Scalar m_thickness;
    Scalar m_stiffness;
    Scalar m_hardness;
    Scalar m_friction;

    bool m_colliding ;
    StaticCollision m_collision ;
    StaticCollision m_collisionAtCreation ;
    mutable bool m_strandP_col_initialized ;
    mutable bool m_strandQ_col_initialized ;

    mutable Vec3x m_p0, m_p1, m_q0, m_q1;
};

struct StrandStrandForceBaseComparator
{
    StrandStrandForceBaseComparator( const ElasticStrand* strand1, int edge1,
            const ElasticStrand* strand2, int edge2 )
    {
        if ( strand1 == strand2 )
        {
            if ( edge1 < edge2 )
            {
                m_edgeP = edge1;
                m_edgeQ = edge2;
            }
            else
            {
                m_edgeP = edge2;
                m_edgeQ = edge1;
            }
        }
        else if ( strand1 > strand2 )
        {
            m_strandP = strand1;
            m_strandQ = strand2;
            m_edgeP = edge1;
            m_edgeQ = edge2;
        }
        else
        {
            m_strandP = strand2;
            m_strandQ = strand1;
            m_edgeP = edge2;
            m_edgeQ = edge1;
        }
    }

    virtual ~StrandStrandForceBaseComparator()
    {
    }

    virtual bool operator()( const std::tr1::shared_ptr<StrandStrandForce>& autre ) = 0;

    const ElasticStrand* m_strandP;
    const ElasticStrand* m_strandQ;
    int m_edgeP;
    int m_edgeQ;
};

struct StrandStrandForceExactComparator: public StrandStrandForceBaseComparator
{
    StrandStrandForceExactComparator( const ElasticStrand* strand1, int edge1,
            const ElasticStrand* strand2, int edge2 ) :
            StrandStrandForceBaseComparator( strand1, edge1, strand2, edge2 )
    {
    }

    virtual ~StrandStrandForceExactComparator()
    {
    }

    bool operator()( const std::tr1::shared_ptr<StrandStrandForce>& other )
    {
        if ( other.use_count() )
            return m_strandP == other->m_strandP && m_edgeP == other->m_edgeP
                    && m_strandQ == other->m_strandQ && m_edgeQ == other->m_edgeQ;
        else
            return false;
    }
};

struct StrandStrandForceProximityComparator: public StrandStrandForceBaseComparator
{
    StrandStrandForceProximityComparator( const ElasticStrand* strand1, int edge1,
            const ElasticStrand* strand2, int edge2 ) :
            StrandStrandForceBaseComparator( strand1, edge1, strand2, edge2 )
    {
    }

    virtual ~StrandStrandForceProximityComparator()
    {
    }

    bool operator()( const std::tr1::shared_ptr<StrandStrandForce>& other )
    {
        if ( other.use_count() )
            return m_strandP == other->m_strandP && abs( m_edgeP - other->m_edgeP ) <= 1
                    && m_strandQ == other->m_strandQ && abs( m_edgeQ - other->m_edgeQ ) <= 1;
        else
            return false;
    }
};

inline bool SegmentsProjOnEachOther( const Vec3x& p0, const Vec3x& p1, const Vec3x& q0,
        const Vec3x& q1, Scalar &s, Scalar &t )
{
    Vec3x dp = p1 - p0; // Direction vector of segment S1
    Vec3x dq = q1 - q0; // Direction vector of segment S2
    Vec3x r = p0 - q0;

    double a = dp.dot( dp ); // Squared length of segment S1, always nonnegative
    double e = dq.dot( dq ); // Squared length of segment S2, always nonnegative
    double f = dq.dot( r );

    double c = dp.dot( r );

    double b = dp.dot( dq );
    double denom = a * e - b * b;

    if ( denom == 0 )
        return false;

    s = ( b * f - c * e ) / denom;
    t = ( b * s + f ) / e;

    static const double margin = 0.1;

    return ( s >= -margin && s <= 1 + margin && t >= -margin && t <= 1 + margin );
}

inline std::pair<int, int> SegmentsSituation( const Vec3x& p0, const Vec3x& p1, const Vec3x& q0,
        const Vec3x& q1 )
{
    Vec3x dp = p1 - p0; // Direction vector of segment S1
    Vec3x dq = q1 - q0; // Direction vector of segment S2
    Vec3x r = p0 - q0;

    double a = dp.dot( dp ); // Squared length of segment S1, always nonnegative
    double e = dq.dot( dq ); // Squared length of segment S2, always nonnegative
    double f = dq.dot( r );

    double c = dp.dot( r );

    double b = dp.dot( dq );
    double denom = a * e - b * b;

    if ( denom == 0 )
        return std::make_pair( 0, 0 );

    Scalar s = ( b * f - c * e ) / denom;
    Scalar t = ( b * s + f ) / e;

    std::pair<int, int> situation;

    if ( s < 0 )
        situation.first = -1;
    else if ( s > 0 )
        situation.first = +1;
    else
        situation.first = 0;

    if ( t < 0 )
        situation.second = -1;
    else if ( t > 0 )
        situation.second = +1;
    else
        situation.second = 0;

    return situation;
}

} /* namespace strandsim */
#endif /* STRANDSTRANDFORCE_HH_ */
