#ifndef STRETCHINGFORCE_HH_
#define STRETCHINGFORCE_HH_

#include "ViscousOrNotViscous.hh"
#include "../Core/ElasticStrand.hh"

template<typename ViscousT = NonViscous>
class StretchingForce: public ForceBase
{
public:
    StretchingForce()
    {}

    virtual ~StretchingForce()
    {}

public:
    static const IndexType s_first = 0; // The first index on which this force can apply
    static const IndexType s_last = 1; // The last index (counting from the end)

    typedef Eigen::Matrix<Scalar, 6, 1> LocalForceType;
    typedef Eigen::Matrix<Scalar, 6, 6> LocalJacobianType;
    typedef VecXx ForceVectorType;

    static std::string getName()
    {
        return ViscousT::getName() + "stretching";
    }

    static Scalar localEnergy( const ElasticStrand& strand, StrandState& geometry,
            const IndexType vtx );

    static void computeLocal( LocalForceType& localF, const ElasticStrand& strand,
            StrandState& geometry, const IndexType vtx );

    static void computeLocal( LocalJacobianType& localJ, const ElasticStrand& strand,
            StrandState& geometry, const IndexType vtx );

    static void addInPosition( ForceVectorType& globalForce, const IndexType vtx,
            const LocalForceType& localForce );

    static void addInPosition( JacobianMatrixType& globalJacobian, const IndexType vtx,
            const LocalJacobianType& localJacobian );

    static void accumulateCurrentE( Scalar& energy, ElasticStrand& strand );
    static void accumulateCurrentF( VecXx& force, ElasticStrand& strand );
    static void accumulateCurrentJ( JacobianMatrixType& Jacobian, ElasticStrand& strand );
};

#endif /* STRETCHINGFORCE_HH_ */
