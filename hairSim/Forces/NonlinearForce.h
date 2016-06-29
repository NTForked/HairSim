
#ifndef NON_LINEAR_FORCE_H_
#define NON_LINEAR_FORCE_H_

#include "../../bogus/Interfaces/MecheEigenInterface.hpp"
#include "../../bogus/Core/ExternalForce.hh"

/*
    Call-back that we will be called every few iterations for the nonLinearSolverWithContacts
    friction problem. Updates the linear system based on the current velocities, and update the friction
    problem.
*/
class NonLinearForce : public bogus::ExternalForce
{
public:
    NonLinearForce( ImplicitStepper &stepper, bogus::MecheFrictionProblem& mecheProblem, unsigned objectId )
        : ExternalForce( objectId ),
          m_stepper( stepper ),
          m_problem ( mecheProblem )
    {}

    virtual ~NonLinearForce() 
    {}

    virtual bool compute( const VecXx& velocities, const VecXx& solverForces )
    // DK: this is what is called from updateExternalForces() in problem.cc which in turn is called by solver.hh in the main loop
    {
        m_stepper.m_futureVelocities = velocities;
        
        bool needsUpdate = m_stepper.updateLinearSystem( solverForces );
        if( needsUpdate )
        { // if update occured, need to inform Bogus/problem solver of new system
            JacobianSolver *M = &m_stepper.linearSolver();

            std::cout << "updating linear system " << m_stepper.m_strand.getGlobalIndex() <<  std::endl;

            Eigen::MatrixXd MecheM( M->matrix().rows(), M->matrix().cols() );
            for( int r = 0; r < MecheM.rows(); ++r ){
                for( int c= 0; c < MecheM.cols(); ++c ){
                    MecheM(r,c) = M->matrix()(r,c);
                }
            }

            m_problem.updateObjectLHS( m_objectID, MecheM );
            m_problem.updateObjectRHS( m_objectID, -m_stepper.rhs() );
        }
        return needsUpdate;
    }

private:
    ImplicitStepper& m_stepper;
    bogus::MecheFrictionProblem& m_problem;
};


#endif