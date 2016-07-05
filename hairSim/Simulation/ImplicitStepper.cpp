#include <limits>

#include "ImplicitStepper.h"
#include "../Strand/DOFScriptingController.h"
#include "../Strand/StrandDynamics.h"
#include "SimulationParameters.h"

#include "../Strand/ElasticStrand.h"
#include "../Forces/StretchingForce.hh"
#include "../Forces/ForceAccumulator.hh"
#include "../Collision/TwistEdgeHandler.h"
#include "../Forces/NonLinearForce.h"

/*
TODO: this file/code is a long way from that,
but it used to be SolveLinear first, and then if that fails
try SolveNonLinear as necessary afterwards....

Now it just goes directly to Nonlinear
....would be faster/better if it was back to intended use
*/



/*
    A dynamic step can trigger a succession of failsafes.
    Here is a summary of the process

    - startStep()

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
    - If update( true ) returns false, which mean the strand is stretching, the StramdImplicitManager
      will try other failsafe to re-solve the contacts and constraints

    - finalize()

  */

ImplicitStepper::ImplicitStepper( ElasticStrand& strand, const SimulationParameters &params ):
    m_futureVelocities( VecXx::Zero( strand.getCurrentDegreesOfFreedom().rows() ) ), 
    m_params( params ),
    m_notSPD( false ), 
    m_usedNonlinearSolver( false ),
    m_strand( strand ), //
    m_dt( 0. ),
    m_nonlinearCallback( NULL ),
    m_newtonIter( 0 )
{
    const Scalar slenderness = m_strand.getTotalRestLength() / 
            ( m_strand.getRadius( 0 ) * m_strand.getRadius( m_strand.getNumVertices() - 1 ) );

    m_inextensibilityThreshold = std::pow( 10., -m_params.m_inextensibility_threshold );
    m_stretchingFailureThreshold = std::pow( 10., -m_params.m_stretching_threshold );
    m_costretchResidualFailureThreshold = std::pow( 10., -m_params.m_costretch_residual_threshold );

    const Scalar alpha = -0.5 * std::log( slenderness ) ;
    m_stretchDamping = std::exp( m_params.m_stretchDamping * m_params.m_stretchDamping * alpha );
}

ImplicitStepper::~ImplicitStepper()
{}

void ImplicitStepper::startStep( const Scalar& dt )
{
    this->m_dt = dt;
    m_notSPD = false;
    m_strand.requireExactJacobian( true );
    m_usedNonlinearSolver = false;
    m_strand.dynamics().computeViscousForceCoefficients( m_dt ); // Maybe overkill to do this each time but at least we can change the stepper's time step without worrying about it.
}

void ImplicitStepper::solveUnconstrained( bool useNonLinearSolver, const bool& penaltyBefore )
{
    if( !penaltyBefore ){
        prepareDynamics();
    }

    if( useNonLinearSolver )
    {
        solveNonLinear();
    }
    else
    {
        solveLinear();

        m_notSPD = m_linearSolver.notSPD();
        if( m_notSPD )
        {
            if ( m_strand.requiresExactJacobian() ) // We are here for the first time, exact Jacobian gives non SPD LHS. Let's try with approximate Jacobian
            {
                std::cerr << "Strand " << m_strand.getGlobalIndex() << " has non spd lhs, using approximate Jacobian" << std::endl;
                m_strand.requireExactJacobian( false );
                solveUnconstrained();
                return;
            }
        }
        m_linearSolver.solve( m_futureVelocities, rhs() );
    }
} // Should call update so see computed changes from m_futureVelocities

void ImplicitStepper::prepareDynamics()
{ // reset so they match start of timestep
    m_strand.setFutureDegreesOfFreedom( m_strand.getCurrentDegreesOfFreedom() );
    m_futureVelocities = m_strand.dynamics().getCurrentVelocities();
}

void ImplicitStepper::solveNonLinear()
{
    StrandDynamics& dynamics = m_strand.dynamics() ;

    Scalar minErr = 1.e99, prevErr = 1.e99;
    m_newtonIter = 0;

    JacobianMatrixType bestLHS;
    VecXx bestRhs, prevRhs;

    Scalar alpha = 0.5;          // Current step length
    const Scalar minAlpha = 0.1; // Minimum step length

    bool foundOneSPD = false;
    m_strand.requireExactJacobian( false );

    // Newton loop -- try to zero-out m_rhs
    for( m_newtonIter = 0; m_newtonIter < m_params.m_maxNewtonIterations; ++m_newtonIter )
    {     
        dynamics.setDisplacements( m_dt * m_futureVelocities );

        prevRhs = m_rhs;
        computeRHS();

        if( m_newtonIter )
        {
            VecXx residual = m_rhs;

            dynamics.getScriptingController()->fixRHS( residual );
            const Scalar err = residual.squaredNorm() / residual.size();

            if( err < minErr || ( !foundOneSPD && !m_linearSolver.notSPD() ) )
            {
                foundOneSPD = !m_linearSolver.notSPD();
                minErr = err;

                if( isSmall( err ) || ( m_newtonIter > 3 && minErr < 1.e-6 ) )
                {
                    m_rhs = prevRhs;
                    break;
                }

                bestLHS = Lhs();
                bestRhs = prevRhs;
            }

            // Decrease or increase the step length based on current convergence
            if( err < prevErr ){
                alpha = std::min( 1.0, 1.5 * alpha );
            }
            else{
                alpha = std::max( minAlpha, 0.5 * alpha );
            }

            prevErr = err;
        }

        computeLHS();

        m_rhs = m_rhs * alpha;

        Lhs().multiply( m_rhs, 1., m_futureVelocities );
        dynamics.getScriptingController()->fixLHSAndRHS( Lhs(), m_rhs, m_dt );
        
        m_linearSolver.store( Lhs() );
        m_notSPD = m_linearSolver.notSPD();
        m_linearSolver.solve( m_futureVelocities, m_rhs );
    }

    // If the non-linear solve failed, returns to the the least problematic step
    if( m_newtonIter == m_params.m_maxNewtonIterations )
    {
        m_rhs = bestRhs;
        Lhs() = bestLHS;

        m_linearSolver.store( Lhs() );
        m_linearSolver.solve( m_futureVelocities, rhs() );
        m_notSPD = m_linearSolver.notSPD();
    }

    m_usedNonlinearSolver = true;
}

void ImplicitStepper::solveLinear()
{
    computeRHS();
    computeLHS();

    Lhs().multiply( m_rhs, 1.0, m_futureVelocities );
    m_strand.dynamics().getScriptingController()->fixLHSAndRHS( Lhs(), m_rhs, m_dt );
    m_linearSolver.store( Lhs() );
}

void ImplicitStepper::computeRHS()
{
    StrandDynamics& dynamics = m_strand.dynamics();

    // start of step and unconstrained
    m_rhs = dynamics.getCurrentVelocities() - m_futureVelocities;
    dynamics.multiplyByMassMatrix( m_rhs );

    const Scalar origKs = m_strand.getParameters().getKs();
    m_strand.getParameters().setKs( m_stretchDamping * origKs );
    dynamics.computeFutureForces( true, true );
    m_strand.getParameters().setKs( origKs );

    VecXx forces = m_strand.getFutureTotalForces();

    m_rhs += forces * m_dt;
}

void ImplicitStepper::computeLHS()
{
    StrandDynamics& dynamics = m_strand.dynamics() ;

    const Scalar origKs = m_strand.getParameters().getKs();
    m_strand.getParameters().setKs( m_stretchDamping * origKs );
    dynamics.computeFutureJacobian( true, true );
    m_strand.getParameters().setKs( origKs );

    JacobianMatrixType& J = m_strand.getTotalJacobian(); // LHS = M - h^2 J
    J *= m_dt * m_dt;

    dynamics.addMassMatrixTo( J );
}
    
VecXx& ImplicitStepper::impulse_rhs()
{
    StrandDynamics& dynamics = m_strand.dynamics();
    m_impulseRhs = m_futureVelocities;
    dynamics.multiplyByMassMatrix( m_impulseRhs );
  
    computeLHS();
    dynamics.getScriptingController()->fixLHSAndRHS( Lhs(), m_impulseRhs, m_dt );
    m_linearSolver.store( Lhs() );

    return m_impulseRhs;
}

// m_futureVelocities has been modified and we want to try and accept the step based on that (or failsafes)
// DK: here's where most of the failsafes kick in -- should modify here to best test algorithm
//     note that this is called by StrandImplicitManager::postProcessFrictionProblem() 
//     where the bool returned here is ignored but m_lastStepWasRejected is set which is checked in StrandImplicitManager::solveCollidingGroup() to see if failsafe is needed
bool ImplicitStepper::update( bool afterConstraints )
{
    bool accept = true;

    VecXx displacements = m_futureVelocities * m_dt;
    m_strand.dynamics().getScriptingController()->enforceDisplacements( displacements ); // enforce any scripting that may have been overriden
    m_strand.dynamics().setDisplacements( displacements );

    // DK: failsafes here:
    const Scalar stretchE = getLineicStretch();
    if( stretchE > m_stretchingFailureThreshold )
    {
        //  Don't bother with very small strands, they would create numerical problems anyway
        if( m_strand.getNumVertices() == 3 && m_strand.getTotalRestLength() < .1 )
        {
            std::cerr << "Strand " << m_strand.getGlobalIndex()
                        << " did bad stuff and is small (len = "
                        << m_strand.getTotalRestLength()
                        << "), reverting to unsimulated pos" << std::endl;

            VecXx futureDofs( m_strand.getCurrentDegreesOfFreedom() ) ;
            m_strand.dynamics().getScriptingController()->setToUnsimulatedPositions( futureDofs, m_strand.getCurrentDegreesOfFreedom() );
            m_strand.setFutureDegreesOfFreedom( futureDofs );
            m_strand.filterFutureGeometryByRestLength( 1.e-6 );
        }
        else
        {
            if( afterConstraints )
            {
                std::cout<< "Strand " << m_strand.getGlobalIndex()
                << " is stretching (" << stretchE << ") after constraints\n";
                accept = false;
            }
            else
            { // Contacts have not been solved yet and we're already stretching ? That's bad.
                // Try again with more favorable starting degrees of freedom
                if( m_strand.requiresExactJacobian() )
                {
                    m_strand.requireExactJacobian( false );
                    solveUnconstrained();
                    return update();
                }
                else
                {
                    accept = false ;
                }
            }

            // Really bad stretching, might need length projection
            if( m_params.m_useLengthProjection && stretchE > 1.0 )
            {
                if( !afterConstraints && !m_usedNonlinearSolver && m_params.m_useNonLinearAsFailsafe ){
                    std::cerr << "Strand " << m_strand.getGlobalIndex()
                                << " is stretching (" << stretchE << "), switch to use non-linear solve " << std::endl;
                    solveUnconstrained( true );
                    return update();
                }
                else{
                    filterGeometryLength( false );
                }
            }

        }
    }

    m_lastStepWasRejected = !accept; // flag so we don't have to catch the call to update(), and still know the outcome of the attempt
    m_strand.getFutureState().freeCachedQuantities();
    return accept;
}

Scalar ImplicitStepper::getLineicStretch()
{
    Scalar stretchE = 0;
    ForceAccumulator<StretchingForce<NonViscous> >::accumulateFuture( stretchE, m_strand );
    stretchE /= m_strand.getParameters().getKs() * m_strand.getTotalRestLength();

    return stretchE;
}

void ImplicitStepper::finalize()
{
    m_strand.dynamics().nanFailSafe();
    m_strand.dynamics().acceptFuture();
}

void ImplicitStepper::prepareForExternalSolve()
{
    if( !m_usedNonlinearSolver )
    { // calls to solve everything if hasnt been done already
        resetStep();
        solveUnconstrained( m_params.m_alwaysUseNonLinear );
    }
}

void ImplicitStepper::resetStep()
{ // if called requires calling SolveUnconstrained again before accepting
    m_usedNonlinearSolver = false;
    prepareDynamics();
}

void ImplicitStepper::filterGeometryLength( bool preStep )
{
    if( preStep ){
        m_strand.filterFutureGeometryByRestLength( m_inextensibilityThreshold );
    }
    else{
        m_strand.filterFutureGeometryByRestLength( m_stretchingFailureThreshold, true );
    }
}

/////
// Nonlinear Callback:
/////

void ImplicitStepper::addNonLinearCallback( bogus::MecheFrictionProblem& mecheProblem, unsigned objectId )
{
    if( m_nonlinearCallback ){
        delete m_nonlinearCallback;
    }
    m_nonlinearCallback = new NonLinearForce( *this, mecheProblem, objectId );
    mecheProblem.addExternalForce( m_nonlinearCallback );
}

// Updates the current Lhs and rhs based of the m_newVelocities guess
bool ImplicitStepper::updateLinearSystem( const VecXx solverForces )
{ //DK: this is where the lineic stretch gets checked and a possible update is applied
    StrandDynamics& dynamics = m_strand.dynamics();
    m_strand.requireExactJacobian( false );
    dynamics.setDisplacements( m_dt * m_futureVelocities ); // updated in NonLinearForce from Bogus, pushing updates 
    
    const Scalar stretchE = getLineicStretch();
    bool needUpdate = stretchE > m_stretchingFailureThreshold;
    Scalar residual = 0;
    if( !needUpdate )
    {
        computeRHS();
        dynamics.getScriptingController()->fixRHS( m_rhs );
        m_usedNonlinearSolver = true;

        //Below should be present, turning off because usually causes needUpdate to occur when stretching is small
        residual = ( m_rhs + solverForces ).squaredNorm() / m_rhs.rows();
        // needUpdate = residual > m_costretchResidualFailureThreshold;

    }
        
    if( needUpdate ) // DK: if residual or stretchE too big.
    { // solveLinear causes less stretching

        solveLinear();
        // solveNonLinear();
    }
        // std::cout << needUpdate << " stretchE: " << stretchE << " | " << m_stretchingFailureThreshold << " ||| residual: " << residual << " | " << m_costretchResidualFailureThreshold << std::endl;
    
    return needUpdate; // DK: now returns true if needs update
}
