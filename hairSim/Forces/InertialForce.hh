#ifndef STRANDSIM_INERTIALFORCE_HH
#define STRANDSIM_INERTIALFORCE_HH

#include "ForceBase.hh"
#include "../Core/BandMatrixFwd.hh"

namespace strandsim {

class InertialForce : public ForceBase
{
public:
    static const IndexType s_first = 0; // The first index on which this force can apply
    static const IndexType s_last = 0; // The last index (counting from the end)

    typedef Vec3x LocalForceType;
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

    static void setFrameAccelerations(
            const Vec3x& Omega, const Vec3x& dOmega_dt,
            const Vec3x& accOrigin ) ;

protected:
    InertialForce();
    virtual ~InertialForce() ;

    static Mat3x s_Omega_Cross ;
    static Mat3x s_Omega_Cross_Cross ;
    static Mat3x s_dOmega_dt_Cross ;
    static Vec3x s_accOrigin ;
};

} // namespace strandsim

#endif // STRANDSIM_INERTIALFORCE_HH
