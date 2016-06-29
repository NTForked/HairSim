#ifndef IMPLICIT_STEPPER_HH
#define IMPLICIT_STEPPER_HH

#include "../Math/SymmetricBandMatrixSolver.h"
#include "../../bogus/Interfaces/MecheEigenInterface.hpp"

#include <tr1/memory>
#include <vector>

namespace bogus
{
    class ExternalForce;
}

class ElasticStrand;
class SimulationParameters;

class ImplicitStepper
{
public:
    ImplicitStepper( ElasticStrand& strand, const SimulationParameters &params );

    virtual ~ImplicitStepper();

    //! Starts a new substep
    void startSubstep( const Scalar& dt );

    //! Solves the unconstrained dynamics, using either linear of non-linear solver
    void solveUnconstrained( bool useNonLinearSolver = true );

    //! Updates the strand degrees of freedom using the velocities stored in m_newVelocities
    /*! Checks if the strand has a high stretch energy. If this is the case and afterContraints if false,
      calls an appropriate failsafe
      \return whether the new velocities are acceptable ( stretch energy is low enough )
    */
    bool update( bool afterConstraints = false );

    //! Rewind the strand and ensures the linear system / constraints are up to date
    void prepareForExternalSolve();

    //! Accept's the new strand's state
    void finalize();

    //! Reverts the strand's state to the one at the beginning of the timestep
    void resetStep();

    //! Updates the current Lhs and rhs based of the m_newVelocities guess
    /*! \return whether the linear system has been updated */
    bool updateLinearSystem( const VecXx solverForces );

    void addNonLinearCallback( bogus::MecheFrictionProblem& problem, unsigned objectId );

    VecXx& impulse_rhs();

    Scalar getDt() const
    { return m_dt; }

    const ElasticStrand& getStrand() const
    { return m_strand; }

    const JacobianMatrixType& Lhs() const
    { return m_strand.getTotalJacobian(); }

    JacobianMatrixType& Lhs()
    { return m_strand.getTotalJacobian(); }

    const VecXx& rhs() const
    { return m_rhs; }
    
    bool notSPD() const
    { return m_notSPD; }

    bool lastStepWasRejected() const
    { return m_lastStepWasRejected; }

    bool usedNonLinearSolver()
    { return m_usedNonlinearSolver; }

    bool refusesMutualContacts() const
    { return notSPD(); }

    JacobianSolver& linearSolver()
    { return m_linearSolver; }

private:

    void prepareDynamics();

    //! Computes linearized dynamics
    void solveLinear();

    //! Computes non-linear dynamics using a Newton algorithm
    void solveNonLinear();

    //! Computes the right-hand-side of the linear system of linearized dynamics at current guess
    void computeRHS();

    //! Computes the left-hand-side of the linear system of linearized dynamics at current guess
    void computeLHS();

    //! Geometric projection of the strand's future dofs to enfore edges rest lengths
    /*! \param preStep  whether this projection is done before or after dynamics */
    void filterGeometryLength( bool preStep );

    //! Returns the value of the stretch energy divided by the length of the rods times its stiffness
    Scalar getLineicStretch();

//////

    VecXx m_futureVelocities;
    VecXx m_impulseChanges;

    SimulationParameters& m_params; // there should only be one, this is a pointer to Scene's simParams
    Scalar m_inextensibilityThreshold;
    Scalar m_stretchingFailureThreshold;
    Scalar m_costretchResidualFailureThreshold;
    Scalar m_stretchDamping;

    JacobianSolver m_linearSolver;

    bool m_notSPD;
    bool m_usedNonlinearSolver;
    bool m_lastStepWasRejected;

    VecXx m_rhs;
    VecXx m_impulseRhs; // for zeroth-order contact 
    
    ElasticStrand& m_strand;
    Scalar m_dt;

    bogus::ExternalForce *m_nonlinearCallback;
    
    unsigned m_newtonIter;

    friend class NonLinearForce;
};

#endif
