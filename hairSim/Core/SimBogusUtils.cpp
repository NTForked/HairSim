bool StrandImplicitManager::assembleBogusFrictionProblem( CollidingGroup& collisionGroup,
        bogus::MecheFrictionProblem& mecheProblem, std::vector<unsigned> &globalIds,
        std::vector<ProximityCollision*> &colPointers, VecXx& vels, VecXx& impulses,
        VecXu& startDofs, VecXu& nDofs, int& numSubSys )
{

    unsigned nContacts = collisionGroup.second.size();

    std::vector<unsigned> nndofs;
    unsigned nSubsystems = collisionGroup.first.size(); 
    numSubSys = nSubsystems;
    globalIds.reserve( nSubsystems );

    nndofs.resize( nSubsystems );
    nDofs.resize( nSubsystems );
    startDofs.resize( nSubsystems );

    unsigned dofCount = 0;
    unsigned subCount = 0;

    std::vector<unsigned> dofIndices;

    // Computes the number of DOFs per strand and the total number of contacts
    for ( IndicesMap::const_iterator it = collisionGroup.first.begin();
            it != collisionGroup.first.end(); ++it )
    {
        startDofs[ subCount ] = dofCount;
        dofIndices.push_back( dofCount );
        const unsigned sIdx = it->first;
        globalIds.push_back( sIdx );
        nDofs[subCount] =  m_steppers[sIdx]->velocities().rows();
        nndofs[subCount] =  m_steppers[sIdx]->velocities().rows();
        nContacts += m_externalContacts[sIdx].size();
        dofCount += m_steppers[sIdx]->velocities().rows();
        ++subCount;
    }

    if( nContacts == 0 ){
        std::cerr << " No contacts EXITING" << std::endl;
        return false; // needs to be false to break out of solvecollidinggroup if statement
    }
    // cout << "contacts found, attempting" << endl;
    // Prepare (rewind) the steppers
#pragma omp parallel for
    for ( unsigned i = 0; i < globalIds.size(); ++i )
    {
        m_steppers[globalIds[i]]->prepareForExternalSolve();
    }

    colPointers.resize( nContacts );

    std::vector < strandsim::SymmetricBandMatrixSolver<double, 10>* > MassMat;
    MassMat.resize( nSubsystems );

    VecXx forces( dofCount );

    vels.resize( dofCount );
    impulses.resize( nContacts * 3 );

    VecXx mu( nContacts );

    std::vector < Eigen::Matrix< double, 3, 3 > > E; E.resize( nContacts );

    VecXx u_frees ( nContacts * 3 );

    int ObjA[ nContacts ];
    int ObjB[ nContacts ];

    std::vector < SparseRowMatx* > H_0; H_0.resize( nContacts );
    std::vector < SparseRowMatx* > H_1; H_1.resize( nContacts );

    unsigned objectId = 0, collisionId = 0, currDof = 0;
    bool oneNonSPD = false;


    // Setting up objects and adding external constraints
    for ( IndicesMap::const_iterator it = collisionGroup.first.begin(); it != collisionGroup.first.end(); ++it )
    {
        const unsigned sIdx = it->first;
        ImplicitStepper& stepper = *m_steppers[sIdx];

        if ( stepper.notSPD() )
            oneNonSPD = true;

        if ( ! m_params.m_useImpulseMethod )
        {
            strandsim::JacobianSolver *M = &stepper.linearSolver();
            MassMat[ objectId++ ] = M;

            if ( m_params.m_alwaysUseNonLinear )
            {
                forces.segment( currDof, stepper.rhs().size() ) = -stepper.rhs();
                currDof += stepper.rhs().size();
            }
            else
            {
                forces.segment( currDof, stepper.impulse_rhs().size() ) = -stepper.impulse_rhs();
                currDof += stepper.impulse_rhs().size();
            }
        }
        else
        {
            // Pass a "pure" mass matrix - impulse problem
            strandsim::JacobianSolver *M = &stepper.massMatrixLinearSolver();
            MassMat[ objectId++ ] = M;

            forces.segment( currDof, stepper.impulse_rhs().size() ) = -stepper.impulse_rhs();
            currDof += stepper.impulse_rhs().size();
        }

        ProximityCollisions & externalCollisions = m_externalContacts[sIdx];
        for ( unsigned i = 0; i < externalCollisions.size(); ++i )
        {
            ProximityCollision& c = externalCollisions[i];
            const Scalar frictionCoeff = c.mu;

            mu[ collisionId ] = frictionCoeff;
            E [ collisionId ] = c.transformationMatrix;
            u_frees.segment<3>( ( collisionId ) * 3 ) = c.objects.first.freeVel - c.objects.second.freeVel;
    // cout << "external ufrees.dotNormal : " << c.normal.dot( u_frees.segment<3>( ( collisionId ) * 3 )  ) << endl;

            ObjA[ collisionId ] = (int) it->second;
            ObjB[ collisionId ] = -1;
            H_0 [ collisionId ] = c.objects.first.defGrad;
            H_1 [ collisionId ] = NULL; // c.objects.second.defGrad;
            impulses.segment<3>( ( collisionId ) * 3 ) = c.force;
            colPointers[collisionId++] = &c;
        }
    }


    // Setting up mutual constraints
#pragma omp parallel for
    for ( unsigned i = 0; i < collisionGroup.second.size(); ++i )
    {
        ProximityCollision& collision = collisionGroup.second[i];
        const int oId1 = collisionGroup.first.find( collision.objects.first.globalIndex )->second;
        const int oId2 = collisionGroup.first.find( collision.objects.second.globalIndex )->second;
        const Scalar frictionCoeff = oneNonSPD ? 0. : collision.mu;
        mu[ collisionId + i ] = frictionCoeff;
        E [ collisionId + i ] = collision.transformationMatrix;
        u_frees.segment<3>( ( collisionId + i ) * 3 ) = collision.objects.first.freeVel - collision.objects.second.freeVel;
    // cout << "ufrees.dotNormal : " << collision.normal.dot( u_frees.segment<3>( ( collisionId + i ) * 3 ) ) << " f:" << collision.objects.first.freeVel.norm() << " s: " << collision.objects.second.freeVel.norm() <<  endl;
        ObjA[ collisionId + i ] = oId1;
        ObjB[ collisionId + i ] = oId2;
        H_0 [ collisionId + i ] = collision.objects.first.defGrad;
        H_1 [ collisionId + i ] = collision.objects.second.defGrad;
        impulses.segment<3>( ( collisionId + i ) * 3 ) = collision.force;
        colPointers[collisionId + i] = &collision;
    }
    assert( collisionId + collisionGroup.second.size() == nContacts );

/*
    cout << "nSubsystems: " << nSubsystems << endl;
    for( int n = 0; n < nndofs.size(); ++n ){
        cout << "nndofs[" << n << "]: " << nndofs[n] << endl; 
    }
    cout << "forces: " << forces << endl;
    cout << "mu: " << mu << endl;
    cout << "nContacts: " << nContacts << endl;
    
    for( int e = 0; e < E.size(); ++e ){
        cout << "E[" << e << "]: " << E[e] << endl; 
    }

    cout << "u_frees: " << u_frees << endl;
    for( int h = 0; h < nContacts; ++h ){
        cout << "ObjA[" << h << "]: " << ObjA[h] << endl; 
    }    
    for( int h = 0; h < nContacts; ++h ){
        cout << "ObjB[" << h << "]: " << ObjB[h] << endl; 
    }    


    for( int h = 0; h < H_0.size(); ++h ){
        cout << "H_0[" << h << "]: " << *(H_0[h]) << endl; 
    }
    for( int h = 0; h < H_1.size(); ++h ){
        cout << "H_1[" << h << "]: " << *(H_1[h]) << endl; 
    }
*/

    mecheProblem.fromPrimal( 
                    nSubsystems,
                    nndofs,
                    MassMat, 
                    forces,
                    nContacts, 
                    mu, 
                    E, 
                    u_frees,
                    ObjA, 
                    ObjB, 
                    H_0, 
                    H_1,
                    dofIndices );

    return true;
}

void StrandImplicitManager::postProcessBogusFrictionProblem( bool accept, CollidingGroup& collisionGroup,
        const bogus::MecheFrictionProblem& , const std::vector<unsigned> &globalIds,
        const std::vector<ProximityCollision*> &colPointers, VecXx& vels, VecXx& impulses,
        VecXu& startDofs, VecXu& nDofs  )
{
    if ( accept )
    {
        assert( globalIds.size() == (unsigned) startDofs.size() );
        assert( globalIds.size() == (unsigned) nDofs.size() );

        m_globalIds.clear();
        m_globalIds = globalIds;

#pragma omp parallel for
        for ( unsigned i = 0; i < globalIds.size(); ++i )
        {
            const unsigned sIdx = globalIds[i];
            const unsigned subSystem = collisionGroup.first.find( sIdx )->second;
// cout << "postProcessBogusFrictionProblem, updating: " << sIdx << endl;
            //UPDATING absolute end of timestep velocity <==> vels 
            m_steppers[sIdx]->newVelocities() = vels.segment( startDofs[ subSystem ], nDofs[ subSystem ] );

            // m_collisionDetector->m_proxyHistory->applyImpulses( m_strands[sIdx], m_steppers[sIdx], true );

            m_steppers[sIdx]->update( true );
        }
    }
    else
    {
        cerr << "StrandImplicitManager::postProcessBogusFrictionProblem... unacceptable result... rewinding " << endl;
#pragma omp parallel for
        for ( unsigned i = 0; i < globalIds.size(); ++i )
        {
            const unsigned sIdx = globalIds[i];
            m_steppers[sIdx]->rewind();
        }
    }

}

bool StrandImplicitManager::solveBogusFrictionProblem( bogus::MecheFrictionProblem& mecheProblem,
        const std::vector<unsigned> &globalIds, bool asFailSafe, bool , VecXx& vels, VecXx& impulses, int& numSubSys )
{

    if ( nonLinearCallbackBogus )
    {
        for ( unsigned i = 0; i < globalIds.size(); ++i )
        {
            m_steppers[globalIds[i]]->addNonLinearCallback( mecheProblem, i );
        }
    }

    bool notGood = true;
    try
    {
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

        // if ( residual > 1e-5 ){
        //     cout << "Bogus residual: " << residual << endl;
        // }
        int numContacts = impulses.size() / 3;
        notGood = residual > std::sqrt( m_params.m_gaussSeidelTolerance ); // This is completely arbitrary
        if ( notGood ){
            cerr << "GS did not converge [ err=" << residual << ", n=" << numContacts << " ] ";
            notGood = false;
        }
    } 
    catch ( std::exception& ex ){
        std::cerr << "StrandImplicitManager::solveBogusFrictionProblem failure... exception" << std::endl;
#ifndef NDEBUG
        throw ex;
#endif
    }

    return !notGood;
}