#include "StrandDynamics.hh"
#include "DOFScriptingController.hh"

#include "../Core/ElasticStrand.hh"
#include "../Core/ElasticStrandUtils.hh"

#include "../Forces/ViscousOrNotViscous.hh"
#include "../Forces/ForceAccumulator.hh"

#include "../Forces/StretchingForce.hh"
#include "../Forces/TwistingForce.hh"
#include "../Forces/BendingForce.hh"
#include "../Forces/GravitationForce.hh"
#include "../Forces/AirDragForce.hh"
#include "../Forces/InertialForce.hh"

#include "../Render/StrandRenderer.hh"
#include "../Utils/TextLog.hh"

StrandDynamics::StrandDynamics( ElasticStrand &strand ) :
        m_strand( strand ),
        m_scriptingController( NULL ), //
        m_futureJacobianUpToDate( false ), //
        m_futureForcesUpToDate( false ), //
        m_DOFmassesUpToDate( false ),
{}

StrandDynamics::~StrandDynamics()
{}

void StrandDynamics::resizeSelf()
{
    const unsigned ndofs = m_strand.getCurrentDegreesOfFreedom().rows();
    m_DOFmasses.resize( ndofs );
}

void StrandDynamics::computeViscousForceCoefficients(Scalar dt)
{
    computeDOFMasses();
    m_strand.getParameters().computeViscousForceCoefficients( dt );
}

void StrandDynamics::computeDOFMasses()
{
    if( m_DOFmassesUpToDate ) return ;

    for( IndexType vtx = 0 ; vtx < m_strand.m_numVertices ; ++vtx )
    {
        m_DOFmasses[4 * vtx + 0] = m_DOFmasses[4 * vtx + 1] = m_DOFmasses[4 * vtx + 2] =
                m_strand.m_vertexMasses[vtx];

        if ( vtx < m_strand.m_numEdges ){
            m_DOFmasses[4 * vtx + 3] = m_strand.getEdgeInertia( vtx );
        }
    }

    m_DOFmassesUpToDate = true ;
}

////////////////////////////////////////////////////////////////////////////////
// Dynamic methods, using viscous forces
////////////////////////////////////////////////////////////////////////////////

void StrandDynamics::computeFutureJacobian( bool withViscous, bool butOnlyForBendingModes )
{
    if ( m_futureJacobianUpToDate )
    {
        return;
    }

    StrandState& futureState = *m_strand.m_futureState ;

    JacobianMatrixType& futureJ = *( futureState.m_totalJacobian );
    futureJ.setZero();

    m_strand.accumulateJ< StretchingForce<NonViscous> > ( futureState ) ;
    m_strand.accumulateJ< TwistingForce<NonViscous> > ( futureState ) ;
    m_strand.accumulateJ< BendingForce<NonViscous> > ( futureState ) ;

    if ( withViscous )
    {
        if ( !butOnlyForBendingModes )
        {
            m_strand.accumulateJ< StretchingForce<Viscous> > ( futureState ) ;
        }
        m_strand.accumulateJ< TwistingForce<Viscous> > ( futureState ) ;
        m_strand.accumulateJ< BendingForce<Viscous> > ( futureState ) ;
        m_strand.accumulateJ< AirDragForce > ( futureState ) ;
    }

    m_strand.accumulateJ< InertialForce > ( futureState ) ;

    futureJ *= -1.0; // To match BASim's sign conventions

    m_futureJacobianUpToDate = true;

    // Free some memory
    futureState.m_hessTwists.free();
    futureState.m_hessKappas.free();
}

void StrandDynamics::computeLHS( Scalar dt, bool withViscous )
{
    computeFutureJacobian( withViscous );
    JacobianMatrixType& LHS = m_strand.getTotalJacobian();
    LHS *= dt * dt;
    addMassMatrixTo( LHS );
    getScriptingController()->fixLHS( LHS ); // Enforce scripted vertices

}

void StrandDynamics::computeFutureForces( bool withViscous, bool butOnlyForBendingModes )
{
    if ( m_futureForcesUpToDate )
    {
        return;
    }

    StrandState& futureState = *m_strand.m_futureState ;

    futureState.m_totalEnergy = 0.0; // NB energy is not going to be used
    VecXx& futureF = futureState.m_totalForce;
    futureF.setZero();

    m_strand.accumulateF< StretchingForce<NonViscous> > ( futureState ) ;
    m_strand.accumulateF< TwistingForce<NonViscous> > ( futureState ) ;
    m_strand.accumulateF< BendingForce<NonViscous> > ( futureState ) ;

    if ( withViscous )
    {
        if ( !butOnlyForBendingModes )
        {
            m_strand.accumulateF< StretchingForce<Viscous> > ( futureState ) ;
        }
        m_strand.accumulateF< TwistingForce<Viscous> > ( futureState ) ;
        m_strand.accumulateF< BendingForce<Viscous> > ( futureState ) ;
        m_strand.accumulateF< AirDragForce > ( futureState ) ;
    }

    m_strand.accumulateF< InertialForce > ( futureState ) ;
    m_strand.accumulateF< GravitationForce > ( futureState ) ;

    m_futureForcesUpToDate = true;
}

void StrandDynamics::computeFutureConservativeEnergy()
{
    StrandState& futureState = *m_strand.m_futureState ;

    m_strand.accumulateE< StretchingForce<NonViscous> > ( futureState ) ;
    m_strand.accumulateE< TwistingForce<NonViscous> > ( futureState ) ;
    m_strand.accumulateE< BendingForce<NonViscous> > ( futureState ) ;
    m_strand.accumulateE< GravitationForce > ( futureState ) ;
}

void StrandDynamics::addMassMatrixTo( JacobianMatrixType& J ) const
{
    for ( int i = 0; i < m_DOFmasses.size(); i++ )
    {
        J( i, i ) += m_DOFmasses[i];
    }
}

const VecXx& StrandDynamics::getDOFMasses() const
{
    return m_DOFmasses;
}

void StrandDynamics::multiplyByMassMatrix( VecXx& tmp ) const
{
    tmp.array() *= m_DOFmasses.array();
}

void StrandDynamics::acceptFuture()
{
    // future will no longer be valid, and current will be set correctly
    m_currentVelocities = m_strand.getFutureDegreesOfFreedom() - m_strand.getCurrentDegreesOfFreedom();
    m_strand.swapStates();

    // m_strand.getFutureDegreesOfFreedom() = m_strand.getCurrentDegreesOfFreedom();
}

void StrandDynamics::nanFailSafe()
{

    std::cerr << "this needs to get changed to operate on future DoFs " << std::endl;
    if( containsNans( m_strand.getCurrentDegreesOfFreedom() ) )
    {
        std::cerr << "Elastic strand " << m_strand.m_globalIndex << " was contaminated by NaNs: reverting to rigid motion" << std::endl;
        m_strand.swapStates();

        const VecXx& currentDOFs = m_strand.getCurrentDegreesOfFreedom();
        VecXx futureDOFs;
        m_scriptingController->computeRigidBodyMotion( futureDOFs, currentDOFs );

        if( containsNans( futureDOFs ) )
        {
            std::cerr << "Elastic strand " << m_strand.m_globalIndex << " has still NaNs in it: reverting to unsimulated vertices" << std::endl;
            m_scriptingController->setToUnsimulatedPositions( futureDOFs, currentDOFs, true );

            if( containsNans( futureDOFs ) )
            {
                std::cerr << "Damn it. " << std::endl;
            }
        }

        m_strand.setFutureDegreesOfFreedom( futureDOFs );

        acceptGuess(); // Re-compute displacements as well, as they were problably NaNised
    }
}
