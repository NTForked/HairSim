#ifndef ELASTICPARAMETERS_HH_
#define ELASTICPARAMETERS_HH_

#include "DependencyNode.hh"

/**
 * Unit: cm
 */
class PhysicalRadius: public DependencyNode< Scalar >
{
public:
    PhysicalRadius( Scalar radius ) :
            DependencyNode< Scalar >( radius )
    {
#ifdef VERBOSE_DEPENDENCY_NODE
        std::cout << "Creating " << name() << ' ' << this << '\n';
#endif

        setClean();
    }

    virtual const char* name() const
    {
        return "PhysicalRadius";
    }

    friend class boost::serialization::access;
    template<class Archive>
    void save( Archive & ar, const unsigned int version ) const
    {
        // Can't use get() or compute() here as boost require save() to be const.
        // So we'll just issue a warning if someone tries to save a dirty value.
        if ( isDirty() )
        {
            std::cerr << "Saving dirty value for " << name() << std::endl;
        }
        ar & m_value;
    }
    template<class Archive>
    void load( Archive & ar, const unsigned int version )
    {
        Scalar value;
        ar & value;
        set( value );
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()

protected    :
    virtual void compute()
    {}
};

/**
 * Unit: no dimension
 */
class BaseRotation: public DependencyNode<Scalar>
{
public:
    BaseRotation( Scalar baseRotation ) :
            DependencyNode<Scalar>( baseRotation )
    {
#ifdef VERBOSE_DEPENDENCY_NODE
        std::cout << "Creating " << name() << ' ' << this << '\n';
#endif

        setClean();
    }

    virtual const char* name() const
    {
        return "BaseRotation";
    }

protected:
    virtual void compute()
    {}
};

/**
 * \brief This contains the bending matrix base, must be multiplied by the appropriate viscous or non-viscous
 * coefficient (with optional interpolation factor).
 *
 * Unit: cm^4
 */
class BendingMatrixBase: public DependencyNode<Mat2x>
{
public:
    BendingMatrixBase( PhysicalRadius& rad, BaseRotation& baseRotation ) :
            DependencyNode<Mat2x>( Mat2x() ), //
            m_physicalRadius( rad ), //
            m_baseRotation( baseRotation )
    {
#ifdef VERBOSE_DEPENDENCY_NODE
        std::cout << "Creating " << name() << ' ' << this << '\n';
#endif

        m_physicalRadius.addDependent( this );
        m_baseRotation.addDependent( this );
    }

    virtual const char* name() const
    {
        return "BendingMatrixBase";
    }

protected:
    virtual void compute()
    {
        const Scalar& radius = m_physicalRadius.get();
        const Scalar baseRotation = m_baseRotation.get();

        Mat2x& B = m_value;
        B( 0, 0 ) = M_PI_4 * radius * cube( radius );
        B( 1, 1 ) = M_PI_4 * radius * cube( radius );
        // rotate cross section by a constant angle
        const Mat2x& rot = Eigen::Rotation2D<Scalar>( baseRotation ).toRotationMatrix();
        B = rot * B * rot.transpose();
        B( 0, 1 ) = B( 1, 0 ) = 0.5 * ( B( 0, 1 ) + B( 1, 0 ) ); // For perfect numerical symmetry

        setDependentsDirty();
    }

    PhysicalRadius& m_physicalRadius;
    BaseRotation& m_baseRotation;
};

/**
 * Unit: dPa = g cm^-1 s^-2
 */
class YoungsModulus: public DependencyNode<Scalar>
{
public:
    YoungsModulus( Scalar youngsModulus ) :
            DependencyNode<Scalar>( youngsModulus )
    {
        setClean();
    }

    virtual const char* name() const
    {
        return "YoungsModulus";
    }

protected:
    virtual void compute()
    {}
};

/**
 * Unit: dPa = g cm^-1 s^-2
 */
class ShearModulus: public DependencyNode<Scalar>
{
public:
    ShearModulus( Scalar shearModulus ) :
            DependencyNode<Scalar>( shearModulus )
    {
        setClean();
    }

    virtual const char* name() const
    {
        return "ShearModulus";
    }

protected:
    virtual void compute()
    {}
};

/**
 * Unit: 10^-5 N = g cm s^-2
 */
class ElasticKs: public DependencyNode<Scalar>
{
public:
    ElasticKs( PhysicalRadius& rad, YoungsModulus& ym ) :
            DependencyNode<Scalar>( std::numeric_limits<Scalar>::signaling_NaN() ), //
            m_physicalRadius( rad ), //
            m_youngsModulus( ym ) //
    {
        m_physicalRadius.addDependent( this );
        m_youngsModulus.addDependent( this );
    }

    virtual const char* name() const
    {
        return "ElasticKs";
    }

protected:
    virtual void compute()
    {
        const Scalar& radius = m_physicalRadius.get();
        const Scalar youngsModulus = m_youngsModulus.get();

        m_value = M_PI * radius * radius * youngsModulus;

        setDependentsDirty();
    }

    PhysicalRadius& m_physicalRadius;
    YoungsModulus& m_youngsModulus;
};

/**
 * Unit: 10^-5 cm^2 N = g cm^3 s^-2
 */
class ElasticKt: public DependencyNode<Scalar>
{
public:
    ElasticKt( PhysicalRadius& rad, ShearModulus& sm ) :
            DependencyNode<Scalar>( std::numeric_limits<Scalar>::signaling_NaN() ), //
            m_physicalRadius( rad ), //
            m_shearModulus( sm )
    {
        m_physicalRadius.addDependent( this );
        m_shearModulus.addDependent( this );
    }

    virtual const char* name() const
    {
        return "ElasticKt";
    }

protected:
    virtual void compute()
    {
        const Scalar& radius = m_physicalRadius.get();
        const Scalar shearModulus = m_shearModulus.get();

        m_value = M_PI_4 * radius * radius
                * ( radius * radius + radius * radius ) * shearModulus;

        setDependentsDirty();
    }

    PhysicalRadius& m_physicalRadius;
    ShearModulus& m_shearModulus;
};

#endif /* ELASTICPARAMETERS_HH_ */
