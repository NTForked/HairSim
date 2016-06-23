#include "StrandImplicitManager.hh"
#include "ImplicitStepper.hh"
#include "StrandDynamicTraits.hh"
#include "../Collision/CollisionDetector.hh"
#include "../Collision/ElementProxy.hh"
#include "../Collision/EdgeFaceIntersection.hh"
#include "../Collision/ContinuousTimeCollision.hh"
#include "../Collision/VertexFaceCollision.hh"
#include "../Collision/EdgeFaceCollision.hh"
#include "../Collision/EdgeEdgeCollision.hh"
#include "../Collision/CollisionUtils.hh"
#include "../Core/ElasticStrand.hh"
#include "../Forces/LevelSetForce.hh"
#include "../Static/ReverseStaticStepper.hh"
#include "../Render/StrandRenderer.hh"
#include "../Utils/SpatialHashMap.hh"
#include "../Utils/LoggingTimer.hh"
#include "../Forces/LevelSetForce.hh"
#include "../Render/Color.hh"

#include "Config.hh"
#include "../../bogus/Interfaces/MecheEigenInterface.hpp"

#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <vector>
#include <fstream>
#include <map>
#include <omp.h>
#include <iostream>

#define SECOND_EDGE_MIN_CONTACT_ABSCISSA 0.0001
#define ALMOST_PARALLEL_COS 0.96592582628 // cos( Pi/12 )

using namespace std;

StrandImplicitManager::StrandImplicitManager( 
        const std::vector<ElasticStrand*>& strands,
        Scalar& dt, 
        const SimulationParameters& params )
: m_dt( dt )
, m_params( params )
, m_strands( strands )
, m_steppers()
, m_collisionDetector( NULL )
, m_hashMap( NULL )
, m_num_ct_hair_hair_col ( 0 )
, m_unconstrained_NewtonItrs( 0 )
{
    accumulateProxies( originalProxies );
    m_externalContacts.resize( m_strands.size() );
    m_collisionDetector = new CollisionDetector( originalProxies, radius );
}

StrandImplicitManager::~StrandImplicitManager()
{
    for( auto stepper = m_steppers.begin(); stepper != m_steppers.end(); ++stepper ){
        delete *stepper;
    }

    for( auto elem = m_elementProxies.begin(); elem != m_elementProxies.end(); ++elem ){
        delete *elem;
    }

    delete m_hashMap;
    m_hashMap = NULL;

    delete m_collisionDetector;
}

void StrandImplicitManager::execute( int substepsInDt )
{
    for( int substep = 0; substep < substepsInDt; ++substep )
    {
        step( substep, m_dt );
    }
}

int hIter;
bool traversalsOn = true;
void StrandImplicitManager::step( unsigned subStepId, Scalar dt )
{

    hIter = 0;
    int hMaxIter = 5;
    bool collisionResolution = true;

    step_prepare( dt );
    m_mutualContacts.clear();
    setupProximityCollisions( dt );

    step_dynamics( subStepId, dt );
    // std::cout  << "cannot access positions until after step_dynamics, otherwise may be reading in wrong q vector" << std::endl;

    m_collisionDetector->m_proxyHistory->m_frozenScene = false;
    m_collisionDetector->m_proxyHistory->trackTunneling = false;

    setupContinuousTimeCollisions(); // should do a first pass where we use regular oldschool collision resolution 
    doContinuousTimeDetection( dt );

    if( collisionResolution ){
        m_collidingGroups.clear();
        m_collidingGroupsIdx.assign( m_strands.size(), -1 );
        step_processCollisions( dt ); // this takes care of current collisions
        step_solveCollisions(); // This is where collisions get solved and strands are finalized() (dV & dX accepted)            
        m_mutualContacts.clear();
    }

    if( hLoop )
    {
        cout << "entering hloop" << endl;
        while( hIter < hMaxIter && !isCollisionInvariantCT( dt ) )
        {
            ++hIter;
            cout << "found unresolved contacts after loop iteration: " << hIter << endl;
            m_collidingGroups.clear();
            m_collidingGroupsIdx.assign( m_strands.size(), -1 );
            step_processCollisions( dt ); // this takes care of current collisions
            step_solveCollisions(); // This is where collisions get solved and strands are finalized() (dV & dX accepted)
            m_mutualContacts.clear();
        }
        cout << "hIter: " << hIter << endl;
    }

    if( trackGeometricRelations )
    {
        m_collisionDetector->clear();
        m_collisionDetector->m_proxyHistory->trackTunneling = true;
        setupContinuousTimeCollisions(); // create and detect loop for missed collisions
        m_collisionDetector->m_proxyHistory->m_frozenScene = true;
        
        if( traversalsOn) {
            ++(m_collisionDetector->m_proxyHistory->m_frozenCheck);
            traversalCheck();  // traverse (mesh optimization)
            ++(m_collisionDetector->m_proxyHistory->m_frozenCheck);    
        }
        std::cout << "avgnumbands: " << m_collisionDetector->m_proxyHistory->tunnelingBands.size() / m_strands.size() << std::endl;
        deleteInvertedProxies(); // delete valid triangles
    }

    m_collisionDetector->clear();
    m_mutualContacts.clear();
    m_time += dt;
}

void StrandImplicitManager::step_prepare( Scalar dt )
{
    for ( unsigned i = 0; i < m_externalContacts.size(); ++i )
    {
        m_externalContacts[i].clear();
    }
}

void StrandImplicitManager::step_dynamics( unsigned subStepId, Scalar substepDt )
{
    // Dynamics system assembly
#pragma omp parallel for schedule(dynamic, 10)
    for ( std::vector<ElasticStrand*>::size_type i = 0; i < m_strands.size(); i++ )
    {
        m_steppers[i]->m_dt = substepDt; // required for checkpointing

        if( penaltyOnce || !trackGeometricRelations || penaltyAfter ) m_collisionDetector->m_proxyHistory->applyImpulses( m_strands[i], m_steppers[i], false ); // false, just setting up m
        else if( !penaltyAfter ) m_collisionDetector->m_proxyHistory->applyImpulses( m_strands[i], m_steppers[i], true );

        m_steppers[i]->startSubstep( subStepId, substepDt );

        m_steppers[i]->solveUnconstrained();
        m_steppers[i]->update();
    }
}

static const unsigned maxObjForOuterParallelism = 8 * omp_get_max_threads();
void StrandImplicitManager::step_processCollisions( Scalar dt )
{
    computeCollidingGroups( m_mutualContacts, dt );

    // Deformation gradients at constraints
    unsigned nExternalContacts = 0;
#pragma omp parallel for reduction ( + : nExternalContacts )
    for ( std::vector<ElasticStrand*>::size_type i = 0; i < m_strands.size(); i++ )
    {
        if ( m_params.m_useDeterministicSolver ){
            std::sort( m_externalContacts[i].begin(), m_externalContacts[i].end() );
        }
        for ( unsigned k = 0; k < m_externalContacts[i].size(); ++k )
        {
            setupDeformationBasis( m_externalContacts[i][k] );
        }

        nExternalContacts += m_externalContacts[i].size();
    }

#pragma omp parallel for
    for ( std::vector<CollidingGroup>::size_type i = 0; i < m_collidingGroups.size(); i++ )
    {
        CollidingGroup& cg = m_collidingGroups[i];

        if ( cg.first.size() < maxObjForOuterParallelism ){
            for ( unsigned k = 0; k < cg.second.size(); ++k )
            {
                setupDeformationBasis( cg.second[k] );
            }
        }
    }

    for ( std::vector<CollidingGroup>::size_type i = 0; i < m_collidingGroups.size(); i++ )
    {
        CollidingGroup& cg = m_collidingGroups[i];

        if ( cg.first.size() >= maxObjForOuterParallelism )
        {
#pragma omp parallel for
            for ( unsigned k = 0; k < cg.second.size(); ++k )
            {
                setupDeformationBasis( cg.second[k] );
            }
        }
    }
}

// DK: solve all collisions and then finalize (relax theta's and accept state update).
void StrandImplicitManager::step_solveCollisions()
{
    // Contact Dynamics solve
#pragma omp parallel for
    for ( std::vector<CollidingGroup>::size_type i = 0; i < m_collidingGroups.size(); ++i )
    {
        if ( m_collidingGroups[i].first.size() <= maxObjForOuterParallelism ){
            solveCollidingGroup( m_collidingGroups[i], false, m_params.m_alwaysUseNonLinear );
        }
    }

    for ( unsigned i = 0; i < m_collidingGroups.size(); ++i )
    {
        if ( m_collidingGroups[i].first.size() > maxObjForOuterParallelism ){
            solveCollidingGroup( m_collidingGroups[i], false, m_params.m_alwaysUseNonLinear );
        }
    }

#pragma omp parallel for
    for ( std::vector<ElasticStrand*>::size_type i = 0; i < m_strands.size(); i++ )
    {
        if ( m_collidingGroupsIdx[i] == -1 )
        {
            if ( needsExternalSolve( i ) ){
                solveSingleObject( i, false, m_params.m_alwaysUseNonLinear );
            }
            else if ( m_params.m_alwaysUseNonLinear && !m_steppers[i]->usedNonLinearSolver() )
            {
                std::cout <<" this gets called" << std::endl;
                m_steppers[i]->rewind();
                m_steppers[i]->solveUnconstrained( true );
                m_steppers[i]->update();
            }
        }
        m_steppers[i]->finalize();
    }
}

void StrandImplicitManager::updateParameters( const SimulationParameters& params )
{
    // Enforce desired or maximum number of threads
    {
        const int numThreads = params.m_numberOfThreads > 0 ? params.m_numberOfThreads : sysconf( _SC_NPROCESSORS_ONLN );
        omp_set_num_threads( numThreads );
    }
    m_params = params;
}

