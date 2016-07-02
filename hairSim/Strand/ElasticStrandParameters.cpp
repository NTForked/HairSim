#include "ElasticStrandParameters.h"
#include "ElasticStrandUtils.h"

ElasticStrandParameters::ElasticStrandParameters() :
        // Dummy default values ElasticStrand's method do not generate FPEs
        // when the parameters are not set yet
        m_density( 0. ), //
        m_viscosity( 0. ), //
        m_airDrag( 0. ), //
        m_rootRadiusMultiplier( 1. ), //
        m_tipRadiusMultiplier( 1. ), //
        m_physicalRadius( 0. ), //
        m_baseRotation( 0. ), //
        m_bendingMatrixBase( m_physicalRadius, m_baseRotation ), //
        m_youngsModulus( 0. ), //
        m_shearModulus( 0. ), //
        m_ks( m_physicalRadius, m_youngsModulus ), //
        m_kt( m_physicalRadius, m_shearModulus )
{}

ElasticStrandParameters::ElasticStrandParameters( const ElasticStrandParameters& other ) :
        m_density( other.m_density ), //
        m_viscosity( other.m_viscosity ), //
        m_airDrag( other.m_airDrag ), //
        m_rootRadiusMultiplier( other.m_rootRadiusMultiplier ), //
        m_tipRadiusMultiplier( other.m_tipRadiusMultiplier ), //
        m_physicalRadius( other.m_physicalRadius.get() ), //
        m_baseRotation( other.m_baseRotation.get() ), //
        m_bendingMatrixBase( m_physicalRadius, m_baseRotation ), //
        m_youngsModulus( other.m_youngsModulus.get() ), //
        m_shearModulus( other.m_shearModulus.get() ), //
        m_ks( m_physicalRadius, m_youngsModulus ), //
        m_kt( m_physicalRadius, m_shearModulus )
{}

ElasticStrandParameters::ElasticStrandParameters( 
    Scalar radiusA,
    Scalar YoungsModulus, 
    Scalar shearModulus, 
    Scalar density, 
    Scalar viscosity, 
    Scalar airDrag,
    Scalar baseRotation ) :
        m_density( density ), //
        m_viscosity( viscosity ), //
        m_airDrag( airDrag ), //
        m_rootRadiusMultiplier( 1. ), //
        m_tipRadiusMultiplier( 1. ), //
        m_physicalRadius( radiusA ), //
        m_baseRotation( baseRotation ), //
        m_bendingMatrixBase( m_physicalRadius, m_baseRotation ), //
        m_youngsModulus( YoungsModulus ), //
        m_shearModulus( shearModulus ), //
        m_ks( m_physicalRadius, m_youngsModulus ), //
        m_kt( m_physicalRadius, m_shearModulus )
{}

ElasticStrandParameters& ElasticStrandParameters::operator=( const ElasticStrandParameters& other )
{
    m_density = other.m_density;
    m_viscosity = other.m_viscosity;
    m_airDrag = other.m_airDrag;
    m_rootRadiusMultiplier = other.m_rootRadiusMultiplier;
    m_tipRadiusMultiplier = other.m_tipRadiusMultiplier;

    // Only the base dependencyNodes need to be copied
    m_physicalRadius.set( other.m_physicalRadius.get() );
    m_baseRotation.set( other.m_baseRotation.get() );
    m_youngsModulus.set( other.m_youngsModulus.get() );
    m_shearModulus.set( other.m_shearModulus.get() );

    return *this;
}

void ElasticStrandParameters::setBaseRotation( Scalar baseRotation )
{
    m_baseRotation.set( baseRotation );
}

void ElasticStrandParameters::setDensity( Scalar density )
{
    m_density = density;
}

void ElasticStrandParameters::setRadius( Scalar radiusA )
{
    m_physicalRadius.set( radiusA );
}

void ElasticStrandParameters::setShearModulus( Scalar shearModulus )
{
    m_shearModulus.set( shearModulus );
}

void ElasticStrandParameters::setYoungsModulus( Scalar youngsModulus )
{
    m_youngsModulus.set( youngsModulus );
}

void ElasticStrandParameters::setViscosity( Scalar viscosity )
{
    m_viscosity = viscosity;
}

void ElasticStrandParameters::setAirDrag( Scalar airDrag )
{
    m_airDrag = airDrag;
}

void ElasticStrandParameters::computeViscousForceCoefficients( Scalar dt )
{
    if( dt == m_dt ) return;

    const Scalar m_radiusA = m_physicalRadius.get();

    // Force coefficients are computed without the varying radius multiplier;
    // correct interpolation will be applied when they are accessed
    m_dt = dt;
    m_viscousKs = M_PI * m_radiusA * m_radiusA * 3 * m_viscosity / dt;
    m_viscousKt = M_PI_4 * m_radiusA * m_radiusA * ( m_radiusA * m_radiusA + m_radiusA * m_radiusA ) * m_viscosity / dt;
    m_viscousBendingCoefficientBase = 3 * m_viscosity / dt;
}

Scalar ElasticStrandParameters::interpolatedRadiusMultiplier( int vtx ) const
{
    const Scalar s = vtx / ( m_numVertices - 1. );
    return ( s * m_tipRadiusMultiplier + ( 1. - s ) * m_rootRadiusMultiplier );
}
