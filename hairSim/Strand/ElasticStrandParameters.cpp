#include "ElasticStrandParameters.hh"
#include "ElasticStrandUtils.hh"

ElasticStrandParameters::ElasticStrandParameters() :
        // Dummy default values ElasticStrand's method do not generate FPEs
        // when the parameters are not set yet
        m_density( 0. ), //
        m_viscosity( 0. ), //
        m_airDrag( 0. ), //
        m_rootRadiusMultiplier( 1. ), //
        m_tipRadiusMultiplier( 1. ), //
        m_ellipticalRadii( 0., 0. ), //
        m_baseRotation( 0. ), //
        m_bendingMatrixBase( m_ellipticalRadii, m_baseRotation ), //
        m_youngsModulus( 0. ), //
        m_shearModulus( 0. ), //
        m_ks( m_ellipticalRadii, m_youngsModulus ), //
        m_kt( m_ellipticalRadii, m_shearModulus )
{}

ElasticStrandParameters::ElasticStrandParameters( const ElasticStrandParameters& other ) :
        m_density( other.m_density ), //
        m_viscosity( other.m_viscosity ), //
        m_airDrag( other.m_airDrag ), //
        m_rootRadiusMultiplier( other.m_rootRadiusMultiplier ), //
        m_tipRadiusMultiplier( other.m_tipRadiusMultiplier ), //
        m_ellipticalRadii( other.m_ellipticalRadii.get().first, other.m_ellipticalRadii.get().second ), //
        m_baseRotation( other.m_baseRotation.get() ), //
        m_bendingMatrixBase( m_ellipticalRadii, m_baseRotation ), //
        m_youngsModulus( other.m_youngsModulus.get() ), //
        m_shearModulus( other.m_shearModulus.get() ), //
        m_ks( m_ellipticalRadii, m_youngsModulus ), //
        m_kt( m_ellipticalRadii, m_shearModulus )
{}

ElasticStrandParameters::ElasticStrandParameters( Scalar radiusA, Scalar radiusB,
        Scalar YoungsModulus, Scalar shearModulus, Scalar density, Scalar viscosity, Scalar airDrag,
        Scalar baseRotation ) :
        m_density( density ), //
        m_viscosity( viscosity ), //
        m_airDrag( airDrag ), //
        m_rootRadiusMultiplier( 1. ), //
        m_tipRadiusMultiplier( 1. ), //
        m_ellipticalRadii( radiusA, radiusB ), //
        m_baseRotation( baseRotation ), //
        m_bendingMatrixBase( m_ellipticalRadii, m_baseRotation ), //
        m_youngsModulus( YoungsModulus ), //
        m_shearModulus( shearModulus ), //
        m_ks( m_ellipticalRadii, m_youngsModulus ), //
        m_kt( m_ellipticalRadii, m_shearModulus )
{
    if( radiusA != radiusB ){
        std::cerr << "RADIUS A does not equal RADIUS B, elliptical strands!" << std::endl;
    }
}

ElasticStrandParameters& ElasticStrandParameters::operator=( const ElasticStrandParameters& other )
{
    m_density = other.m_density;
    m_viscosity = other.m_viscosity;
    m_airDrag = other.m_airDrag;
    m_rootRadiusMultiplier = other.m_rootRadiusMultiplier;
    m_tipRadiusMultiplier = other.m_tipRadiusMultiplier;

    // Only the base dependencyNodes need to be copied
    m_ellipticalRadii.set( other.m_ellipticalRadii.get() );
    m_baseRotation.set( other.m_baseRotation.get() );
    m_youngsModulus.set( other.m_youngsModulus.get() );
    m_shearModulus.set( other.m_shearModulus.get() );

    m_ks( m_ellipticalRadii, m_youngsModulus );
    m_kt( m_ellipticalRadii, m_shearModulus );

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

void ElasticStrandParameters::setRadii( Scalar radiusA, Scalar radiusB )
{
    m_ellipticalRadii.set( std::make_pair( radiusA, radiusB ) );
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

    const std::pair<Scalar, Scalar>& radii = m_ellipticalRadii.get();
    const Scalar m_radiusA = radii.first;
    const Scalar m_radiusB = radii.second;

    // Force coefficients are computed without the varying radius multiplier;
    // correct interpolation will be applied when they are accessed
    m_dt = dt;
    m_viscousKs = M_PI * m_radiusA * m_radiusB * 3 * m_viscosity / dt;
    m_viscousKt = M_PI_4 * m_radiusA * m_radiusB * ( m_radiusA * m_radiusA + m_radiusB * m_radiusB ) * m_viscosity / dt;
    m_viscousBendingCoefficientBase = 3 * m_viscosity / dt;
}

Scalar ElasticStrandParameters::interpolatedRadiusMultiplier( int vtx ) const
{
    const Scalar s = vtx / ( m_numVertices - 1. );
    return ( s * m_tipRadiusMultiplier + ( 1. - s ) * m_rootRadiusMultiplier );
}

// Returns the multiplier e( theta ) such that the ellipse equation is
// r( theta ) = e( theta ) * majorRadius
// majorRadius being radiusA
Scalar ElasticStrandParameters::getRadiusShrinkAtAngle( const Scalar angle ) const
{
    const std::pair<Scalar, Scalar>& radii = m_ellipticalRadii.get();
    const Scalar m_radiusA = radii.first;
    const Scalar m_radiusB = radii.second;

    const Scalar a2 = m_radiusA * m_radiusA;
    const Scalar b2 = m_radiusB * m_radiusB;
    const Scalar ct = std::cos( angle );
    return m_radiusA / std::sqrt( b2 + ( a2 - b2 ) * ct * ct );

    return radiusShrinkAtAngle( m_radiusA / m_radiusB, angle );
}

Scalar ElasticStrandParameters::radiusShrinkAtAngle( const Scalar ratio, const Scalar angle )
{
    const Scalar ct = std::cos( angle );
    return ratio / std::sqrt( 1. + ( ratio * ratio - 1. ) * ct * ct );
}
