#include "Simulation.h"
#include "../Collision/CollisionDetector.h"
#include "../Collision/CollisionUtils/SpatialHashMap.hh"
#include <omp.h>

using namespace std;

Simulation::Simulation( 
        const std::vector<ElasticStrand*>& strands,
        SimulationParameters& params,
        const std::vector< TriMesh* >& meshes )
: m_collisionDetector( NULL )
, m_params( params )
, m_strands( strands )
, m_steppers()
, m_hashMap( NULL )
{
    std::vector< ElementProxy* > originalProxies;
    accumulateProxies( originalProxies, meshes );
    m_externalContacts.resize( m_strands.size() );
    m_collisionDetector = new CollisionDetector( originalProxies );
}

Simulation::~Simulation()
{
    for( auto stepper = m_steppers.begin(); stepper != m_steppers.end(); ++stepper ){
        delete *stepper;
    }

    delete m_hashMap;
    m_hashMap = NULL;

    delete m_collisionDetector;
}

int hIter;
bool hLoop = false;
bool traversalsOn = false;
bool trackGeometricRelations = true;
bool penaltyAfter = true;
void Simulation::step( const Scalar& dt )
{
    hIter = 0;
    int hMaxIter = 5;
    bool collisionResolution = true;

    step_prepare( dt );
    step_dynamics( dt );

    if( collisionResolution ){
        gatherProximityRodRodCollisions( dt );
        detectContinuousTimeCollisions(); // should do a first pass where we use regular oldschool collision resolution 
        preProcessContinuousTimeCollisions( dt );

        step_processCollisions( dt ); // this takes care of current collisions
        step_solveCollisions(); // This is where collisions get solved and strands are finalized() (dV & dX accepted)            
    }

    if( hLoop )
    {
        cout << "entering hloop" << endl;
        while( hIter < hMaxIter && !isCollisionInvariantCT( dt ) )
        {
            ++hIter;
            cout << "found unresolved contacts after loop iteration: " << hIter << endl;
            step_processCollisions( dt ); // this takes care of current collisions
            step_solveCollisions(); // This is where collisions get solved and strands are finalized() (dV & dX accepted)
        }
        cout << "hIter: " << hIter << endl;
    }

    // No further dynamics on strands beyond this point, just book-keeping

    if( trackGeometricRelations )
    {
        m_collisionDetector->clear();
        m_collisionDetector->m_proxyHistory->trackTunneling = true;
        detectContinuousTimeCollisions(); // create and detect loop for missed collisions
        m_collisionDetector->m_proxyHistory->m_frozenScene = true;
        
        deleteInvertedProxies( penaltyAfter, !trackGeometricRelations );
    }

    step_finish();
}

void Simulation::step_prepare( Scalar dt )
{
    for ( unsigned i = 0; i < m_externalContacts.size(); ++i )
    {
        m_externalContacts[i].clear();
    }
    m_mutualContacts.clear();

    m_collisionDetector->m_proxyHistory->m_frozenScene = false;
    m_collisionDetector->m_proxyHistory->trackTunneling = false;
}

void Simulation::step_dynamics( Scalar dt )
{
    // Dynamics system assembly
#pragma omp parallel for schedule(dynamic, 10)
    for( std::vector< ElasticStrand* >::size_type i = 0; i < m_strands.size(); ++i )
    {
        m_steppers[i]->setDt( dt ); // required for checkpointing, this needs to be here so long as anything occurs before startSubstep
        m_collisionDetector->m_proxyHistory->applyImpulses( m_strands[i], m_steppers[i], !penaltyAfter );

        m_steppers[i]->startStep( dt );
        m_steppers[i]->solveUnconstrained( true, !penaltyAfter );
        m_steppers[i]->update();
    }
}

static const unsigned maxObjForOuterParallelism = 8 * omp_get_max_threads();
void Simulation::step_processCollisions( Scalar dt )
{
    m_collidingGroups.clear();
    m_collidingGroupsIdx.assign( m_strands.size(), -1 );
    computeCollidingGroups( m_mutualContacts );

    // Deformation gradients at constraints
    unsigned nExternalContacts = 0;
#pragma omp parallel for reduction ( + : nExternalContacts )
    for ( std::vector<ElasticStrand*>::size_type i = 0; i < m_strands.size(); i++ )
    {
        std::sort( m_externalContacts[i].begin(), m_externalContacts[i].end() );

        for( unsigned k = 0; k < m_externalContacts[i].size(); ++k )
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

    for( std::vector<CollidingGroup>::size_type i = 0; i < m_collidingGroups.size(); i++ )
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

void Simulation::step_solveCollisions()
{
    // Contact solve
#pragma omp parallel for
    for( std::vector<CollidingGroup>::size_type i = 0; i < m_collidingGroups.size(); ++i )
    {
        if( m_collidingGroups[i].first.size() <= maxObjForOuterParallelism ){
            solveCollidingGroup( m_collidingGroups[i], false, m_params.m_alwaysUseNonLinear );
        }
    }

    for( unsigned i = 0; i < m_collidingGroups.size(); ++i )
    {
        if( m_collidingGroups[i].first.size() > maxObjForOuterParallelism ){
            solveCollidingGroup( m_collidingGroups[i], false, m_params.m_alwaysUseNonLinear );
        }
    }

#pragma omp parallel for
    for( std::vector<ElasticStrand*>::size_type i = 0; i < m_strands.size(); ++i )
    {
        if ( m_collidingGroupsIdx[i] == -1 ){
            if ( needsExternalSolve( i ) ){
                solveOnlyStrandExternal( i, false, m_params.m_alwaysUseNonLinear );
            }
            else if( m_params.m_alwaysUseNonLinear && !m_steppers[i]->usedNonLinearSolver() ){
                m_steppers[i]->resetStep();
                m_steppers[i]->solveUnconstrained( true );
                m_steppers[i]->update();
            }
        }
    }
    m_mutualContacts.clear();    
}

void Simulation::step_finish()
{
#pragma omp parallel for
    for( std::vector<ElasticStrand*>::size_type i = 0; i < m_strands.size(); ++i )
    {
        m_steppers[i]->finalize(); // Accept and finish with strand motion
    }
    m_collisionDetector->clear();
    m_mutualContacts.clear();
}

void Simulation::updateParameters( const SimulationParameters& params )
{
    { // Enforce desired or maximum number of threads
        const int numThreads = params.m_numberOfThreads > 0 ? params.m_numberOfThreads : sysconf( _SC_NPROCESSORS_ONLN );
        omp_set_num_threads( numThreads );
    }
    m_params = params;
}

