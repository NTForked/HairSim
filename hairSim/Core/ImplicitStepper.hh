#ifndef STRANDSIM_IMPLICITSTEPPER_HH
#define STRANDSIM_IMPLICITSTEPPER_HH

#include "../Core/StepperBase.hh"
#include "../Utils/SymmetricBandMatrixSolver.hh"
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

    //! Solves the unconstrained dynamics, using either linear of non-linear solver
    void solveUnconstrained( bool useNonLinearSolver = true ); //false );
    //! Updates the strand degrees of freedom using the velocities stored in m_newVelocities
    /*! Checks if the strand has a high stretch energy. If this is the case and afterContraints if false,
      calls an appropriate failsafe
      \return whether the new velocities are acceptable ( stretch energy is low enough )
    */
    bool update(bool afterConstraints = false );

    //! Rewind the strand and ensures the linear system / constraints are up to date
    void prepareForExternalSolve() ;

    //! Accept's the new strand's state
    void finalize();

    //! Scales the dynamics' linear sytem by s
    void scale( const Scalar s );

    //! Reverts the strand's state to the one at the beginning of the timestep
    void rewind();

    Scalar getDt() const
    { return m_dt; }

    //! Starts a new substep
    void startSubstep( unsigned id, Scalar m_dt );

    const ElasticStrand& getStrand() const
    {
        return m_strand;
    }

    const JacobianMatrixType& ImplicitStepper::Lhs() const
    {
        return m_strand.getTotalJacobian();
    }

    JacobianMatrixType& ImplicitStepper::Lhs()
    {
        return m_strand.getTotalJacobian();
    }

    const VecXx& rhs() const
    {
        return m_rhs;
    }
    
    VecXx& impulse_rhs() ;
    

    /// H velocities, like displacements, should not exist... make everything depend on strandstate
    VecXx & velocities()
    {
        return m_velocities;
    }
    VecXx & newVelocities()
    {
        return m_newVelocities;
    }

    bool notSPD() const
    {
        return m_notSPD;
    }

    bool lastStepWasRejected() const
    {
        return m_lastStepWasRejected ;
    }

    bool usedNonLinearSolver()
    {
        return m_usedNonlinearSolver ;
    }

    bool refusesMutualContacts() const
    {
        return notSPD();
    }

    JacobianSolver& linearSolver ()
    {
        return m_linearSolver;
    }
    
    JacobianSolver& massMatrixLinearSolver ()
    {
        return m_massMatrix_linearSolver;
    }
    
    //    JacobianSolver& complianceLinearSolver ()
    //    {
    //        return m_compliance_linearSolver;
    //    }

    //! Updates the current Lhs and rhs based of the m_newVelocities guess
    /*! \return whether the linear system has been updated */
    bool updateLinearSystem( const VecXx solverForces ) ;

    void addNonLinearCallback( bogus::MecheFrictionProblem& problem, unsigned objectId ) ;

    unsigned numIters()
    {
        return m_newtonIter;
    }
    
    JacobianMatrixType m_massMatrix; // really just mass. N.B. YacFS overloads mass matrix to stepper's Jacobian
    VecXx m_velocities;
    VecXx m_impulseChanges;

private:
    //! Intialized future degrees of freedom, setup frames, optionally initialize length constraints
    void prepareDynamics() ;

    //! Computes linearized dynamics
    void solveLinear();
    //! Computes non-linear dynamics using a Newton algorithm
    void solveNonLinear();

    //! Computes the right-hand-side of the linear system of linearized dynamics at current guess
    void computeRHS() ;
    //! Computes the left-hand-side of the linear system of linearized dynamics at current guess
    void computeLHS() ;

    //! Copy reference frames from current state to future state
    void setupFuturesFrames() ;

    //! Geometric projection of the strand's future dofs to enfore edges rest lengths
    /*! \param preStep  whether this projection is done before or after dynamics */
    void filterGeometryLength( bool preStep );

    //! Returns the value of the stretch energy divided  by the length of the rods times its stiffness
    Scalar getLineicStretch();

    SimulationParameters* m_params; // there should only be one, this is a pointer to Scene's simParams
    Scalar m_inextensibilityThreshold;
    Scalar m_stretchingFailureThreshold;
    Scalar m_costretchResidualFailureThreshold;
    Scalar m_stretchDamping;

    JacobianSolver m_linearSolver;
    
    JacobianSolver m_massMatrix_linearSolver; // for zeroth-order contact resolve
    
    //    JacobianSolver m_compliance_linearSolver; 
    //    JacobianMatrixType m_complianceMatrix; 

    bool m_notSPD;
    bool m_usedNonlinearSolver;
    bool m_linearSystemIsDirty;
    bool m_lastStepWasRejected;

    VecXx m_newVelocities;
    VecXx m_rhs;
    VecXx m_impulseRhs; // for zeroth-order contact 
    VecXx m_projectionDisplacements;

public:   
    
    ElasticStrand& m_strand;
    Scalar m_dt;

    bogus::ExternalForce *m_nonlinearCallback ;

private:
    
    unsigned m_newtonIter;

};

#endif // STRANDSIM_IMPLICITSTEPPER_HH
