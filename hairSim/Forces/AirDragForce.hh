#ifndef AIRDRAGFORCE_HH
#define AIRDRAGFORCE_HH

#include "ForceBase.hh"
#include "../Math/BandMatrixFwd.h"
#include "../Strand/ElasticStrand.h"

class AirDragForce : public ForceBase
{
public:
    static const IndexType s_first = 0; // The first index on which this force can apply
    static const IndexType s_last = 0; // The last index (counting from the end)

    typedef Vec3 LocalForceType;
    typedef Mat3x LocalJacobianType;
    typedef VecXx ForceVectorType;

    virtual std::string getName() const
    {
        return "air drag" ;
    }

    static void computeLocal( LocalForceType& localF, const ElasticStrand& strand,
            const StrandState& geometry, const IndexType vtx );

    static void computeLocal( LocalJacobianType& localJ, const ElasticStrand& strand,
            const StrandState& geometry, const IndexType vtx );

    static void addInPosition( ForceVectorType& globalForce, const IndexType vtx,
            const LocalForceType& localForce );

    static void addInPosition( JacobianMatrixType& globalJacobian, const IndexType vtx,
            const LocalJacobianType& localJacobian );

    static void setFrameVelocities(
            const Vec3& Omega,
            const Vec3& velOrigin ) ;

protected:
    AirDragForce();
    virtual ~AirDragForce() ;

    static Vec3 s_velOrigin ;
    static Mat3x s_Omega_Cross ;
};

#endif
