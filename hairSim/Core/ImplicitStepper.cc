#include <limits>

#include "ImplicitStepper.hh"
#include "DOFScriptingController.hh"
#include "StrandDynamicTraits.hh"
#include "SimulationParameters.hh"

#include "../Core/ElasticStrand.hh"
#include "../Utils/TextLog.hh"
#include "../Core/ElasticStrand.hh"
#include "../Forces/StretchingForce.hh"
#include "../Forces/ForceAccumulator.hh"
#include "../Collision/TwistEdgeHandler.hh"

#include "../../bogus/Interfaces/MecheEigenInterface.hpp"
#include "../../bogus/Core/ExternalForce.hh"

/*
    A dynamic step can trigger a succession of failsafes.
    Here is a summary of the process

    - startSubstep()

    Before contacts / hard constraints : solveUnconstrained() and update( false )
    - First try with linearized dynamics and exact jacobian.
    - If the jacobian is not SPD or the resulting step makes the strand stretch,
      use the non-exact jacobian, which sould be better conditionned.
    - If the jacobian is still not SPD or the strand is still stretching,
      fall back to the unconstraint non-linear solver ( solveNonLinear() )
    - If the strand is still badly stretching, use a geometric projection

    Optionally for contacts/hard-constraints:
    - Call prepareForExternalSolve()
    - Then either update( true ) to accept the step or rewind() to discard it
    - If update( true ) returns false, which mean the strand is stretching, the  StramdImplicitManager
      will try other failsafe to re-solve the contacts and constraints

    - finalize()

  */

ImplicitStepper::ImplicitStepper( ElasticStrand& strand,
        const SimulationParameters &params ) :
        m_params( params ),
        m_notSPD( false ), 
        m_usedNonlinearSolver( false ),
        m_dt( 0. ),
        m_velocities( VecXx::Zero( strand.getCurrentDegreesOfFreedom().rows() ) ), 
        m_newVelocities( VecXx::Zero( strand.getCurrentDegreesOfFreedom().rows() ) ), 
        m_projectionDisplacements ( VecXx::Zero( m_velocities.rows() ) ),
        m_strand( strand ), //
        m_nonlinearCallback( NULL ),
        m_newtonIter( 0 )
{
    const Scalar slenderness = m_strand.getTotalRestLength()
            / ( m_strand.getRadiusA( 0 ) * m_strand.getRadiusB( m_strand.getNumVertices() - 1 ) );

    m_inextensibilityThreshold = std::pow( 10., -m_params.m_inextensibility_threshold );
    m_stretchingFailureThreshold = std::pow( 10., -m_params.m_stretching_threshold );
    m_costretchResidualFailureThreshold = std::pow( 10., -m_params.m_costretch_residual_threshold );

    const Scalar alpha = - .5 * std::log( slenderness ) ;
    m_stretchDamping = std::exp( m_params.m_stretchDamping * m_params.m_stretchDamping * alpha ) ;
    
    if ( m_params.m_useImpulseMethod )
    {
        // set up mass matrix for zeroth order solve
        const unsigned ndofs = m_strand.getCurrentDegreesOfFreedom().rows();
        m_massMatrix.resize( ndofs, ndofs );
        m_massMatrix.setZero();
    }
    //    else if ( !m_params.m_alwaysUseNonLinear )
    //    {
    //        const unsigned ndofs = m_strand.getCurrentDegreesOfFreedom().rows();
    //        m_complianceMatrix.resize( ndofs, ndofs );
    //        m_complianceMatrix.setZero();
    //    }
}

ImplicitStepper::~ImplicitStepper()
{}


void ImplicitStepper::startSubstep( unsigned id, Scalar dt )
{
    this->m_dt = dt;
    m_notSPD = false;

    m_strand.requireExactJacobian( true );

    m_usedNonlinearSolver = false ;

    m_strand.dynamics().computeViscousForceCoefficients( m_dt ); // Maybe an overkill to do this each time but at least we can change the stepper's time step without worrying about it.
    
    if ( m_params.m_useImpulseMethod )
    {
        // DK: for now init massMatrix here
        // (need to do this after computeViscousForceCoefficients since this where mass is initialized for now -- awkward.)
        m_massMatrix.setZero();
        m_strand.dynamics().addMassMatrixTo( m_massMatrix );
        m_strand.dynamics().getScriptingController()->fixLHS( m_massMatrix );
        m_massMatrix_linearSolver.store( m_massMatrix );
    }

    // Only required if we read from checkpoint at previous t, but no easy way to check that
    m_strand.dynamics().getDisplacements() += (m_impulseChanges * m_dt);
    m_velocities = m_strand.dynamics().getDisplacements() / m_dt;
}

void ImplicitStepper::solveLinear()
{
    StrandDynamicTraits& dynamics = m_strand.dynamics() ;

    computeRHS();
    computeLHS();

    Lhs().multiply( m_rhs, 1., m_newVelocities );
    dynamics.getScriptingController()->fixLHSAndRHS( Lhs(), m_rhs, m_dt );

    m_linearSolver.store( Lhs() );
}

void ImplicitStepper::setupFuturesFrames()
{
    m_strand.getFutureState().m_referenceFrames1.set(
            m_strand.getCurrentState().m_referenceFrames1.get() ); // We need the old ones to compute the new ones
    m_strand.getFutureState().m_referenceFrames2.set(
            m_strand.getCurrentState().m_referenceFrames2.get() ); // We need the old ones to compute the new ones
    m_strand.getFutureState().m_referenceFrames1.getPreviousTangents() =
            m_strand.getCurrentState().m_referenceFrames1.getPreviousTangents(); // Can we avoid the copies? // FIXME LATER
    m_strand.getFutureState().m_referenceTwists.set(
            m_strand.getCurrentState().m_referenceTwists.get() );
}

void ImplicitStepper::computeRHS()
{
    StrandDynamicTraits& dynamics = m_strand.dynamics() ;

            // start of step and unconstrained
    m_rhs = m_velocities - m_newVelocities ;
    dynamics.multiplyByMassMatrix( m_rhs );

    const Scalar origKs = m_strand.getParameters().getKs();
    m_strand.getParameters().setKs( m_stretchDamping * origKs );

    dynamics.computeFutureForces( true, true );

    m_strand.getParameters().setKs( origKs );

    VecXx forces = m_strand.getFutureTotalForces();

    m_rhs += forces * m_dt;
}
    
VecXx& ImplicitStepper::impulse_rhs()
{
    StrandDynamicTraits& dynamics = m_strand.dynamics() ;
    m_impulseRhs = m_newVelocities;
    dynamics.multiplyByMassMatrix( m_impulseRhs );
  
    if ( m_params.m_useImpulseMethod )
    {
      dynamics.getScriptingController()->enforceVelocities( m_newVelocities, m_dt );
    }
    else
    {
      computeLHS();
      dynamics.getScriptingController()->fixLHSAndRHS( Lhs(), m_impulseRhs, m_dt );
      m_linearSolver.store( Lhs() );
    }

    return m_impulseRhs;
}

void ImplicitStepper::computeLHS()
{
    StrandDynamicTraits& dynamics = m_strand.dynamics() ;

    const Scalar origKs = m_strand.getParameters().getKs();
    m_strand.getParameters().setKs( m_stretchDamping * origKs );

    dynamics.computeFutureJacobian( true, true );

    m_strand.getParameters().setKs( origKs );

    JacobianMatrixType& J = m_strand.getTotalJacobian(); // LHS = M - h^2 J
    J *= m_dt * m_dt;

    dynamics.addMassMatrixTo( J );
}

void ImplicitStepper::scale( const Scalar s )
{
    if ( s != 1. )
    {
        m_linearSolver.setScaling( s );
        m_rhs *= s;
    }
}

void ImplicitStepper::prepareDynamics()
{ // sets m_newVelocites to unconstrained quess

    // rest future state and initial guess
    m_strand.setFutureDegreesOfFreedom( m_strand.getCurrentDegreesOfFreedom() );
    
    if ( m_params.m_usePreFilterGeometry )
    {
        filterGeometryLength( true );
        setupFuturesFrames();
    }

    StrandDynamicTraits& dynamics = m_strand.dynamics() ;
    
    if ( m_params.m_usePreFilterGeometry )
        m_newVelocities = m_projectionDisplacements / m_dt;
    else
        m_newVelocities = dynamics.getDisplacements() / m_dt;
    
    dynamics.getDisplacements().setZero();
}

void ImplicitStepper::solveUnconstrained( bool useNonLinearSolver )
{
    prepareDynamics() ;

    if( useNonLinearSolver )
    {
        solveNonLinear();

    } else {

        solveLinear() ;

        m_notSPD = m_linearSolver.notSPD();


        if ( m_notSPD )
        {
            if ( m_strand.requiresExactJacobian() ) // We are here for the first time, exact Jacobian gives non SPD LHS. Let's try with approximate Jacobian
            {
                DebugStream( g_log, "" ) << "Strand " << m_strand.getGlobalIndex()
                                           << " has non spd lhs, using approximate Jacobian";
                m_strand.requireExactJacobian( false );
                solveUnconstrained();
                return ;
            }
            else // We are here for the second time, LHS is definitely not SPD
            {

                CopiousStream( g_log, "" ) << "Strand " << m_strand.getGlobalIndex()
                                          << " has definitely non spd lhs, solving alone without friction";

            }
        }

        m_linearSolver.solve( newVelocities(), rhs() );
    }
}

void ImplicitStepper::solveNonLinear()
{

    StrandDynamicTraits& dynamics = m_strand.dynamics() ;

    Scalar minErr = 1.e99, prevErr = 1.e99 ;
    m_newtonIter = 0 ;

    JacobianMatrixType bestLHS ;
    VecXx bestRhs, prevRhs ;

    Scalar alpha = .5 ;          // Current step length
    const Scalar minAlpha = .1 ; // Minimum step length

    bool foundOneSPD = false ;
    m_strand.requireExactJacobian( false ) ;

    // Newton loop -- try to zero-out m_rhs
    for( m_newtonIter = 0 ; m_newtonIter < m_params.m_maxNewtonIterations ; ++ m_newtonIter )
    {
        dynamics.getDisplacements() = m_dt * m_newVelocities ;
        const VecXx tentativeDofs = m_strand.getCurrentDegreesOfFreedom() + dynamics.getDisplacements() ;
        m_strand.setFutureDegreesOfFreedom( tentativeDofs );

        prevRhs = m_rhs ;
        computeRHS();

        if( m_newtonIter )
        {
            VecXx residual = m_rhs ;

            dynamics.getScriptingController()->fixRHS( residual );
            const Scalar err = residual.squaredNorm() / residual.size() ;

            if( err < minErr || ( !foundOneSPD && !m_linearSolver.notSPD() ) )
            {
                foundOneSPD = !m_linearSolver.notSPD() ;
                minErr = err ;

                if( isSmall( err ) || ( m_newtonIter > 3 && minErr < 1.e-6 ) )
                {
                    m_rhs = prevRhs ;
                    break ;
                }

                bestLHS = Lhs() ;
                bestRhs = prevRhs ;

            }

            // Decrease or increase the step length based on current convergence
            if ( err < prevErr ) {
                alpha = std::min( 1., 1.5*alpha ) ;
            } else {
                alpha = std::max( minAlpha, .5*alpha ) ;
            }

            prevErr = err ;
        }

        computeLHS();


        m_rhs = m_rhs * alpha  ;

        Lhs().multiply( m_rhs, 1., m_newVelocities );
        dynamics.getScriptingController()->fixLHSAndRHS( Lhs(), m_rhs, m_dt );
        
        //        if ( !m_params.m_alwaysUseNonLinear ) //&& m_complianceMatrix.norm() == 0. )
        //        {
        //            m_complianceMatrix = Lhs();
        //            m_compliance_linearSolver.store( m_complianceMatrix );
        //        }
        
        m_linearSolver.store( Lhs() );
        m_notSPD = m_linearSolver.notSPD() ;
        m_linearSolver.solve( newVelocities(), rhs() );


//        std::cout << m_newtonIter << " Residual: " << prevErr  << " SPD " << !m_linearSolver.notSPD()
//                  << " alpha " << alpha << std::endl ;
//        std::cout << residual.transpose() << std::endl ;

    }

    // If the non-linear solve failed, returns to the the least problematic step
    if( m_newtonIter == m_params.m_maxNewtonIterations )
    {
        m_rhs = bestRhs ;
        Lhs() = bestLHS ;

        m_linearSolver.store( Lhs() );
        m_linearSolver.solve( newVelocities(), rhs() );
        m_notSPD = m_linearSolver.notSPD() ;
    }

    m_usedNonlinearSolver = true ;

    DebugStream( g_log, "" ) << "Nonlinear solve for strand " << m_strand.getGlobalIndex() <<
                 " iter " << m_newtonIter << " err " << minErr << " ; SPD: " << !m_notSPD ;

}

    
// Updates the current Lhs and rhs based of the m_newVelocities guess
bool ImplicitStepper::updateLinearSystem( const VecXx solverForces )
//DK: this is where the lineic stretch gets checked and a possible update is applied
{
    StrandDynamicTraits& dynamics = m_strand.dynamics() ;
    
    m_strand.requireExactJacobian( false );
    
    const VecXx tentativeDofs = m_strand.getCurrentDegreesOfFreedom()
    + m_dt * m_newVelocities ;
    m_strand.setFutureDegreesOfFreedom( tentativeDofs );
    
    dynamics.getDisplacements() = m_dt * m_newVelocities ;
    
    const Scalar stretchE = getLineicStretch();
    
    bool needUpdate = stretchE > m_stretchingFailureThreshold ;
    
    ContactStream( g_log, "GS" ) << "Strand " << m_strand.getGlobalIndex()
    << " stretch is : " << stretchE << " ( with max stretch : " << m_stretchingFailureThreshold << " )";
    
    if ( !needUpdate )
    {
        computeRHS();
        dynamics.getScriptingController()->fixRHS( m_rhs );
        
        m_usedNonlinearSolver = true ;
        const Scalar residual = ( m_rhs + solverForces ).squaredNorm() / m_rhs.rows() ;
        
        needUpdate = residual > m_costretchResidualFailureThreshold ;
        
        ContactStream( g_log, "GS" ) << "Strand " << m_strand.getGlobalIndex() << " residual is " << residual  << " ( with max residual : " << m_costretchResidualFailureThreshold << " )";
    }
    
    ContactStream( g_log, "GS" ) << "Strand " << m_strand.getGlobalIndex() << " needUpdate = " << needUpdate;
    
    if( needUpdate ) // DK: if residual or stretchE too big.
    {
        solveLinear();
    }
    
    return needUpdate ; // DK: now returns true if needs update
}
    

    // Call-back that we will be called every few iterations for the nonLinearSolverWithContacts
    // friction problem. Updates the linear system based on the current velocities, and update the friction
    // problem.
class NonLinearForce : public bogus::ExternalForce
{
    public:

    NonLinearForce( ImplicitStepper &stepper, bogus::MecheFrictionProblem& mecheProblem, unsigned objectId )
        : ExternalForce( objectId ),
          m_stepper( stepper ),
          m_problem ( mecheProblem )
    {}

    virtual ~NonLinearForce() {}

    virtual bool compute ( const VecXx& velocities, const VecXx& solverForces )
    // DK: this is what is called from updateExternalForces() in problem.cc which in turn is called by solver.hh in the main loop
    {
        m_stepper.newVelocities() = velocities ;
        
        bool needsUpdate = m_stepper.updateLinearSystem( solverForces );
        if( needsUpdate )
        { // if update occured, need to inform Bogus/problem solver of new system
            strandsim::JacobianSolver *M = &m_stepper.linearSolver();

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
        ImplicitStepper& m_stepper ;
        bogus::MecheFrictionProblem& m_problem ;
} ;

void ImplicitStepper::addNonLinearCallback( bogus::MecheFrictionProblem& mecheProblem, unsigned objectId )
{
    if( m_nonlinearCallback ) delete m_nonlinearCallback;
    m_nonlinearCallback = new NonLinearForce( *this, mecheProblem, objectId ) ;
    mecheProblem.addExternalForce( m_nonlinearCallback ); // DK:
}


Scalar ImplicitStepper::getLineicStretch()
{
    Scalar stretchE = 0;
    ForceAccumulator<StretchingForce<NonViscous> >::accumulateFuture( stretchE, m_strand );
    stretchE /= m_strand.getParameters().getKs() * m_strand.getTotalRestLength();

    return stretchE ;
}

// Updates the strand using m_newVelocities and check for stretching. If stretching, calls
// the appropriate fail-safe or returns false
// DK: here's where most of the failsafes kick in -- should modify here to best test algorithm
//     note that this is called by StrandImplicitManager::postProcessFrictionProblem() 
//     where the bool returned here is ignored but m_lastStepWasRejected is set which is checked in StrandImplicitManager::solveCollidingGroup() to see if failsafe is needed
bool ImplicitStepper::update( bool afterConstraints )
{
    // newVelocities are the absolute final end of timestep velocities used to update the implicit euler step
    VecXx displacements = m_newVelocities * m_dt;
    m_strand.dynamics().getScriptingController()->enforceDisplacements( displacements );
// 
    // end of step positions = start of step position + displacement due to velocity
    m_strand.setFutureDegreesOfFreedom( m_strand.getCurrentDegreesOfFreedom() + displacements );

    bool accept = true ;

    // DK: failsafes here:
    const Scalar stretchE = getLineicStretch();
    if ( stretchE > m_stretchingFailureThreshold  )
    {

        if ( m_params.m_failure_testing_on && stretchE > 0.25 && afterConstraints )
        {
            //InfoStream( g_log, "" )
            std::cout<< "\033[31;1m Strand " << m_strand.getGlobalIndex()
            << " is stretching (" << stretchE << ") after constraints  \033[m ; exiting \n";
            // exit(1);
        }
        
        //  Don't bother with very small strands, they would create numerical problems anyway
        if( m_strand.getNumVertices() == 3 && m_strand.getTotalRestLength() < .1 )
        {
            WarningStream( g_log, "" ) << "Strand " << m_strand.getGlobalIndex()
                                       << " did bad stuff and is small (len="
                                       << m_strand.getTotalRestLength()
                                       << "), reverting to unsimulated pos";

            VecXx futureDofs( m_strand.getCurrentDegreesOfFreedom() ) ;
            m_strand.dynamics().getScriptingController()->setToUnsimulatedPositions( futureDofs, m_strand.getCurrentDegreesOfFreedom() );
            m_strand.setFutureDegreesOfFreedom( futureDofs );
            m_strand.filterFutureGeometryByRestLength( 1.e-6 );

        } else {

            if ( afterConstraints )
            {
                //ContactStream( g_log, "" ) << "\033[31;1m Strand " << m_strand.getGlobalIndex()
                //<< " is stretching (" << stretchE << ") after constraints"; // \033[m";
                InfoStream( g_log, "" ) << "\033[31;1m Strand " << m_strand.getGlobalIndex()
                << " is stretching (" << stretchE << ") after constraints  \033[m";
                std::cout<< "\033[31;1m Strand " << m_strand.getGlobalIndex()
                << " is stretching (" << stretchE << ") after constraints  \033[m ; \n";
                accept = false ;
            }
            else // Contacts have not been solved yet and we're already stretching ? That's bad.
            {
                // Try again with more favorable starting degrees of freedom

                if( m_strand.requiresExactJacobian() )
                {

                    DebugStream( g_log, "" ) << "Strand " << m_strand.getGlobalIndex()
                                               << " is stretching (" << stretchE << "), using approximate Jacobian";
                    m_strand.requireExactJacobian( false );

                    solveUnconstrained();
                    return update() ;

                } else {
                    DebugStream( g_log, "" ) << "Strand " << m_strand.getGlobalIndex()
                                               << " is still stretching (" << stretchE << "), may need non-linear solve " ;
                    accept = false ;
                }

            }

            // DK: (length projection)
            // Really bad stretching, might need geomtric projection
            if ( m_params.m_useLengthProjection && stretchE > 1. )
            {
                if( !afterConstraints && !m_usedNonlinearSolver && m_params.m_useNonLinearAsFailsafe ) {

                    CopiousStream( g_log, "" ) << "Strand " << m_strand.getGlobalIndex()
                                               << " is still stretching (" << stretchE << "), using non-linear solve ";
                    solveUnconstrained( true );
                    return update() ;

                } else {
                    filterGeometryLength( false );
                    CopiousStream( g_log, "" ) << "Projected strand " << m_strand.getGlobalIndex()
                                          << " -> max depl: " << m_projectionDisplacements.lpNorm<Eigen::Infinity>();
                }
            }

        }
    }

    m_lastStepWasRejected = !accept ;

    // update displacements and swap states
    m_strand.dynamics().acceptGuess(); // DK: here's where the solution gets transfered from futurestate to currentstate (note *not* in finalize()) and then current state is possibly accepted as "good"
    // now future is start of step, and current is end of step

    m_strand.getFutureState().freeCachedQuantities();

    return accept;
}

// DK: relax thetas in finalize is yet another failsafe here...
void ImplicitStepper::finalize()
{
    if ( m_params.m_useRelaxTheta )
    {
        // Put back tentative state into futureState
        m_strand.swapStates();

        //        // Temporary projection before relaxing thetas
        //        filterGeometryLength( true );

        // Relax
        try
        {
            m_strand.relaxThetas();
        } catch ( std::exception& ex )
        {
            ErrorStream( g_log, "" ) << "Strand " << m_strand.getGlobalIndex()
                    << " got an exception in relaxThetas(): " << ex.what();
            m_strand.resetThetas();
        }
        m_strand.getFutureState().freeCachedQuantities();

        // Undo projection
        //        m_strand.setFutureDegreesOfFreedom( m_strand.getFutureDegreesOfFreedom() - m_projectionDisplacements );
        
        m_strand.swapStates();
    }

    m_strand.dynamics().nanFailSafe();

    // Finite difference computation of velocities
    newVelocities() = m_strand.dynamics().getDisplacements() / m_dt;
    m_velocities = m_newVelocities;
}

void ImplicitStepper::prepareForExternalSolve()
{
    rewind();

    if( m_usedNonlinearSolver )
    {
        solveUnconstrained( );
    }

    updateConstraints();
}

void ImplicitStepper::rewind()
{
    change this to set future equal to current or something less potentially dangerous
    m_strand.swapStates();
}

void ImplicitStepper::filterGeometryLength( bool preStep )
{
    const VecXx orig = m_strand.getFutureDegreesOfFreedom();

    if ( preStep )
    {
        m_strand.filterFutureGeometryByRestLength( m_inextensibilityThreshold );
    }
    else
    {
        m_strand.filterFutureGeometryByRestLength( m_stretchingFailureThreshold, true );
    }

    m_projectionDisplacements = m_strand.getFutureDegreesOfFreedom() - orig;
}
