#ifndef VISCOUSORNOTVISCOUS_HH_
#define VISCOUSORNOTVISCOUS_HH_

#include "../Strand/ElasticStrand.h"

// These classes are taken as template arguments for the internal forces,
// indicating whether we want the non-viscous or the viscous version.
// The forces call their ViscousT's static methods returning the appropriate
// stiffness and "rest shape" (the actual rest-shape for non-viscous or the
// shape at the beginning of time step for viscous).

class NonViscous
{
protected:
    NonViscous()
    {}

    virtual ~NonViscous()
    {}

public:
    static std::string getName()
    {
        return "";
    }

    static Scalar bendingCoefficient( const ElasticStrand& strand, int vtx )
    {
        return strand.m_parameters.bendingCoefficient( vtx );
    }
    static Mat2x bendingMatrix( const ElasticStrand& strand, int vtx )
    {
        return strand.m_parameters.bendingMatrix( vtx );
    }

    static const Vec2 kappaBar( const ElasticStrand& strand, int vtx )
    {
        return strand.m_restKappas[vtx];
    }

    static Scalar kt( const ElasticStrand& strand, int vtx )
    {
        return strand.m_parameters.getKt( vtx );
    }

    static Scalar thetaBar( const ElasticStrand& strand, int vtx )
    {
        return strand.m_restTwists[vtx];
    }

    static Scalar ks( const ElasticStrand& strand, int vtx )
    {
        return strand.m_parameters.getKs( vtx );
    }

    static Scalar ellBar( const ElasticStrand& strand, int vtx )
    {
        return strand.m_restLengths[vtx];
    }

    class NonDissipativeForce
    {};
};

class Viscous
{
protected:
    Viscous()
    {}

    virtual ~Viscous()
    {}

public:
    static std::string getName()
    {
        return "viscous ";
    }

    static Scalar bendingCoefficient( const ElasticStrand& strand, int vtx )
    {
        return strand.m_parameters.viscousBendingCoefficient( vtx );
    }

    static Mat2x bendingMatrix( const ElasticStrand& strand, int vtx )
    {
        return strand.m_parameters.viscousBendingMatrix( vtx );
    }

    static const Vec2 kappaBar( const ElasticStrand& strand, int vtx )
    {
        return strand.m_currentState->m_kappas[vtx];
    }

    static Scalar kt( const ElasticStrand& strand, int vtx )
    {
        return strand.m_parameters.getViscousKt( vtx );
    }

    static Scalar thetaBar( const ElasticStrand& strand, int vtx )
    {
        return strand.m_currentState->m_twists[vtx];
    }

    static Scalar ks( const ElasticStrand& strand, int vtx )
    {
        return strand.m_parameters.getViscousKs( vtx );
    }

    static Scalar ellBar( const ElasticStrand& strand, int vtx )
    {
        return strand.m_currentState->m_lengths[vtx];
    }

    class DissipativeForce
    {};
};

#endif /* VISCOUSORNOTVISCOUS_HH_ */
