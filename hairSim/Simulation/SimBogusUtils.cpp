#include "Simulation.h"
#include "ImplicitStepper.h"

bool Simulation::assembleBogusFrictionProblem( 
        CollidingGroup& collisionGroup,
        bogus::MecheFrictionProblem& mecheProblem, 
        std::vector<unsigned> &globalIds,
        VecXx& vels, 
        VecXx& impulses,
        VecXu& startDofs, 
        VecXu& nDofs )
{

    unsigned nContacts = collisionGroup.second.size();

    std::vector<unsigned> nndofs;
    unsigned nSubsystems = collisionGroup.first.size(); 
    globalIds.reserve( nSubsystems );

    nndofs.resize( nSubsystems );
    nDofs.resize( nSubsystems );
    startDofs.resize( nSubsystems );

    unsigned dofCount = 0;
    unsigned subCount = 0;
    std::vector<unsigned> dofIndices;

    // Computes the number of DOFs per strand and the total number of contacts
    for( IndicesMap::const_iterator it = collisionGroup.first.begin();
         it != collisionGroup.first.end(); 
         ++it )
    {
        startDofs[ subCount ] = dofCount;
        dofIndices.push_back( dofCount );
        const unsigned sIdx = it->first;
        globalIds.push_back( sIdx );
        nDofs[subCount] = m_steppers[sIdx]->m_futureVelocities.rows();
        nndofs[subCount] = m_steppers[sIdx]->m_futureVelocities.rows();
        nContacts += m_externalContacts[sIdx].size();
        dofCount += m_steppers[sIdx]->m_futureVelocities.rows();
        ++subCount;
    }

    if( nContacts == 0 ){
        return false; // needs to be false to break out of solvecollidinggroup if statement
    }

    // Prepare the steppers
#pragma omp parallel for
    for ( unsigned i = 0; i < globalIds.size(); ++i )
    { // make sure steppers are ready for us to read their info/members... solve if not yet solved, reset where necessary
        m_steppers[ globalIds[i] ]->prepareForExternalSolve();
    }

    std::vector < SymmetricBandMatrixSolver<double, 10>* > MassMat;
    MassMat.resize( nSubsystems );

    VecXx forces( dofCount );

    vels.resize( dofCount );
    impulses.resize( nContacts * 3 );
    impulses.setZero();

    VecXx mu( nContacts );

    std::vector < Eigen::Matrix< double, 3, 3 > > E; E.resize( nContacts );

    VecXx u_frees( nContacts * 3 );
    u_frees.setZero();

    int ObjA[ nContacts ];
    int ObjB[ nContacts ];

    std::vector < SparseRowMatx* > H_0; H_0.resize( nContacts );
    std::vector < SparseRowMatx* > H_1; H_1.resize( nContacts );

    unsigned objectId = 0, collisionId = 0, currDof = 0;

    // Setting up objects and adding external constraints
    for ( IndicesMap::const_iterator it = collisionGroup.first.begin(); it != collisionGroup.first.end(); ++it )
    {
        const unsigned sIdx = it->first;
        ImplicitStepper& stepper = *m_steppers[sIdx];

        JacobianSolver *M = &stepper.linearSolver();
        MassMat[ objectId++ ] = M;

        if( m_params.m_alwaysUseNonLinear )
        {
            forces.segment( currDof, stepper.rhs().size() ) = -stepper.rhs();
            currDof += stepper.rhs().size();
        }
        else
        {
            forces.segment( currDof, stepper.impulse_rhs().size() ) = -stepper.impulse_rhs();
            currDof += stepper.impulse_rhs().size();
        }

        CollidingPairs & externalCollisions = m_externalContacts[sIdx];
        for ( unsigned i = 0; i < externalCollisions.size(); ++i, ++collisionId )
        {
            CollidingPair& c = externalCollisions[i];

            mu[ collisionId ] = c.m_mu;
            E [ collisionId ] = c.m_transformationMatrix;
            // u_frees.segment<3>( ( collisionId ) * 3 ) = c.objects.first.freeVel - c.objects.second.freeVel;
            ObjA[ collisionId ] = (int) it->second;
            ObjB[ collisionId ] = -1;
            H_0 [ collisionId ] = c.objects.first.defGrad;
            H_1 [ collisionId ] = NULL;
        }
    }

    // Setting up RodRod constraints
#pragma omp parallel for
    for ( unsigned i = 0; i < collisionGroup.second.size(); ++i )
    {
        CollidingPair& collision = collisionGroup.second[i];
        const int oId1 = collisionGroup.first.find( collision.objects.first.globalIndex )->second;
        const int oId2 = collisionGroup.first.find( collision.objects.second.globalIndex )->second;

        mu[ collisionId + i ] = collision.m_mu;
        E [ collisionId + i ] = collision.m_transformationMatrix;
        // u_frees.segment<3>( ( collisionId + i ) * 3 ) = collision.objects.first.freeVel - collision.objects.second.freeVel;
        ObjA[ collisionId + i ] = oId1;
        ObjB[ collisionId + i ] = oId2;
        H_0 [ collisionId + i ] = collision.objects.first.defGrad;
        H_1 [ collisionId + i ] = collision.objects.second.defGrad;
    }
    assert( collisionId + collisionGroup.second.size() == nContacts );

    mecheProblem.fromPrimal( 
                    nSubsystems,
                    nndofs,
                    MassMat, 
                    forces,
                    nContacts, 
                    mu, 
                    E, 
                    u_frees, // this is what we want the velocity to be given resting contact, set it to zero?
                    ObjA, 
                    ObjB, 
                    H_0, 
                    H_1,
                    dofIndices );

    return true;
}

bool Simulation::solveBogusFrictionProblem( 
        bogus::MecheFrictionProblem& mecheProblem,
        const std::vector<unsigned> &globalIds, 
        bool asFailSafe, 
        bool nonLinear, 
        VecXx& vels, 
        VecXx& impulses )
{
    if( nonLinear )
    {
        for ( unsigned i = 0; i < globalIds.size(); ++i )
        {
            m_steppers[globalIds[i]]->addNonLinearCallback( mecheProblem, i );
        }
    }

    const double residual = mecheProblem.solve( 
                            impulses,   // impulse guess and returned impulse
                            vels,       // returned velocities
                            0,     // max number of threads, else OpenMP default
                            0.0,   // tolerance
                            0,     // max iterations
                            false, // static problem
                            0.0,   // regularization
                            false, // use infinity norm
                            false, // use projected gradient
                            0 ); // cadoux iters

    bool failed = residual > std::sqrt( m_params.m_gaussSeidelTolerance ); // arbitrary tolerance
    if( failed ){
        std::cerr << "GS did not converge [ err=" << residual << ", numContacts=" << impulses.size() / 3 << " ] " << std::endl;
    }

    return !failed;
}

void Simulation::postProcessBogusFrictionProblem( 
        bool accept, 
        CollidingGroup& collisionGroup,
        const bogus::MecheFrictionProblem& problem, 
        const std::vector<unsigned> &globalIds,
        VecXx& vels, 
        VecXx& impulses,
        VecXu& startDofs, 
        VecXu& nDofs )
{
    if( accept )
    {
        assert( globalIds.size() == (unsigned) startDofs.size() );
        assert( globalIds.size() == (unsigned) nDofs.size() );

        m_globalIds.clear();
        m_globalIds = globalIds;

#pragma omp parallel for
        for( unsigned i = 0; i < globalIds.size(); ++i )
        {
            const unsigned sIdx = globalIds[i];
            const unsigned subSystem = collisionGroup.first.find( sIdx )->second;

            m_steppers[sIdx]->m_futureVelocities = vels.segment( startDofs[ subSystem ], nDofs[ subSystem ] );
            m_steppers[sIdx]->update( true );
        }
    }
}

