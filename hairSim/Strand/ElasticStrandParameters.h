#ifndef ELASTICSTRANDPARAMETERS_HH_
#define ELASTICSTRANDPARAMETERS_HH_

#include "Definitions.hh"
#include "../Dependencies/BendingProducts.hh"

#include <boost/serialization/access.hpp>
#include <boost/serialization/version.hpp>

class ElasticStrandParameters
{
public:
    ElasticStrandParameters();

    ElasticStrandParameters( const ElasticStrandParameters& other );

    ElasticStrandParameters( Scalar radiusA, Scalar YoungsModulus,
            Scalar shearModulus, Scalar density, Scalar viscosity, Scalar airDrag,
            Scalar baseRotation );

    ElasticStrandParameters& operator=( const ElasticStrandParameters& other );

    Scalar getBaseRotation() const
    {
        return m_baseRotation.get();
    }
    Scalar getDensity() const
    {
        return m_density;
    }
    Scalar getKs( int vtx ) const
    {
        const Scalar interpol = interpolatedRadiusMultiplier( vtx );
        return interpol * interpol * m_ks.get();
    }
    Scalar getKt( int vtx ) const
    {
        const Scalar interpol = interpolatedRadiusMultiplier( vtx );
        const Scalar interpol2 = interpol * interpol;

        return interpol2 * interpol2 * m_kt.get();
    }
    Scalar getRadius( int vtx ) const
    {
        return interpolatedRadiusMultiplier( vtx ) * m_physicalRadius.get();
    }

    Scalar getShearModulus() const
    {
        return m_shearModulus.get();
    }
    Scalar getYoungsModulus() const
    {
        return m_youngsModulus.get();
    }
    Scalar getViscosity() const
    {
        return m_viscosity;
    }
    Scalar getAirDrag() const
    {
        return m_airDrag;
    }

    Scalar bendingCoefficient( int vtx ) const
    {
        const Scalar interpol = interpolatedRadiusMultiplier( vtx );
        const Scalar interpol2 = interpol * interpol;

        return interpol2 * interpol2 * m_youngsModulus.get();
    }

    Mat2x bendingMatrix( int vtx ) const
    {
        return bendingCoefficient( vtx ) * m_bendingMatrixBase.get();
    }

    void setBaseRotation( Scalar baseRotation );
    void setDensity( Scalar density );
    void setRadius( Scalar radiusA );
    void setShearModulus( Scalar shearModulus );
    void setYoungsModulus( Scalar YoungsModulus );
    void setViscosity( Scalar viscosity );
    void setAirDrag( Scalar airDrag );

    //!TODO  Since we dont' want to save dirty values, we need a way to compute them before
    // void computeDependencies();

    void computeViscousForceCoefficients( Scalar dt );

    void setRootRadiusMultiplier( double rootRadiusMultiplier )
    {
        m_rootRadiusMultiplier = rootRadiusMultiplier;
    }

    void setTipRadiusMultiplier( double tipRadiusMultiplier )
    {
        m_tipRadiusMultiplier = tipRadiusMultiplier;
    }

    Scalar viscousBendingCoefficient( int vtx ) const
    {
        const Scalar interpol = interpolatedRadiusMultiplier( vtx );
        const Scalar interpol2 = interpol * interpol;

        return interpol2 * interpol2 * m_viscousBendingCoefficientBase;
    }

    Mat2x viscousBendingMatrix( int vtx ) const
    {
        return viscousBendingCoefficient( vtx ) * m_bendingMatrixBase.get();
    }

    Scalar getViscousKs( int vtx ) const
    {
        const Scalar interpol =  interpolatedRadiusMultiplier( vtx );

        return interpol * interpol * m_viscousKs;
    }

    Scalar getViscousKt( int vtx ) const
    {
        const Scalar interpol =  interpolatedRadiusMultiplier( vtx );
        const Scalar interpol2 = interpol * interpol;

        return interpol2 * interpol2 * m_viscousKt;
    }

    Scalar getDt() const
    {
        return m_dt;
    }

    IndexType getNumVertices() const
    {
        return m_numVertices;
    }
    IndexType& getNumVertices()
    {
        return m_numVertices;
    }

    Scalar getKs() const
    {
        return m_ks.get();
    }

    void setKs( const Scalar ks )
    {
        m_ks.set( ks );
    }

    BendingMatrixBase& getBendingMatrixBase()
    {
        return m_bendingMatrixBase;
    }

    const Mat2x& bendingMatrixBase() const
    {
        return m_bendingMatrixBase.get();
    }

    Scalar &dt()
    {
        return m_dt;
    }

private:
    friend class CollisionParameters;
    friend class boost::serialization::access;

    Scalar interpolatedRadiusMultiplier( int vtx ) const;

    template<class Archive> 
    void serialize( Archive & ar, const unsigned int version )
    {
        ar & m_physicalRadius;
        ar & m_youngsModulus;
        ar & m_shearModulus;
        ar & m_density;
        ar & m_baseRotation;

        // Versioning -- increase BOOST_CLASS_VERSION below each time new properties are added
        if ( version > 0 )
        {
            ar & m_viscosity;
            if ( version > 1 )
            {
                ar & m_airDrag;
            }
            if ( version > 2 )
            {
                ar & m_numVertices;
                ar & m_rootRadiusMultiplier;
                ar & m_tipRadiusMultiplier;
                ar & m_bendingMatrixBase;
            }
        }

        ar & m_ks;
        ar & m_kt;
    }

    // Size of the strand
    IndexType m_numVertices;

    double m_density;
    double m_viscosity;
    double m_airDrag;

    double m_rootRadiusMultiplier;
    double m_tipRadiusMultiplier;

    // Computed viscous force coefficients. computeViscousForceCoefficients(dt) MUST be called at the beginning of each time step, in case something has changed.
    double m_viscousBendingCoefficientBase;
    Scalar m_viscousKt;
    Scalar m_viscousKs;

    // Dependencies
    mutable PhysicalRadius m_physicalRadius;
    mutable BaseRotation m_baseRotation;
    mutable BendingMatrixBase m_bendingMatrixBase;
    mutable YoungsModulus m_youngsModulus;
    mutable ShearModulus m_shearModulus;
    mutable ElasticKs m_ks;
    mutable ElasticKt m_kt;

    Scalar m_dt; // This really doesn't belong here, mainly for convinience
};


BOOST_CLASS_VERSION( strandsim::ElasticStrandParameters, 3 )

#endif /* ELASTICSTRANDPARAMETERS_HH_ */
