#include "SimulationUtils.h"

void StrandImplicitManager::updateParameters( const SimulationParameters& params )
{
    // Enforce desired or maximum number of threads
    {
        const int numThreads = params.m_numberOfThreads > 0 ? params.m_numberOfThreads : sysconf( _SC_NPROCESSORS_ONLN );
        omp_set_num_threads( numThreads );
    }
    m_params = params;
}

void SimulationUtils::accumulateProxies( originalProxies )
{
    TwistEdge* last = NULL;
    std::vector< ElementProxy* > originalProxies;
    // Store the stepper for each strand
    // Also store edge proxies for collision detection
    for ( std::vector<ElasticStrand*>::const_iterator strand = m_strands.begin(); strand != m_strands.end(); ++strand )
    {
        double strandRadius = (*strand)->m_collisionRadius < 0 ? radius : (*strand)->m_collisionRadius;

        m_steppers.push_back( new ImplicitStepper( **strand, m_params ) );

        // if( m_steppers.size() == m_strands.size() ) radius /= 1.5;

        m_steppers.back()->m_strand.rod_data->m_renderer->m_strandRadius = strandRadius;

        for ( int vtx = 0; vtx < ( *strand )->getNumEdges(); ++vtx )
        {
            const CollisionParameters& cp = ( *strand )->collisionParameters();

            TwistEdge* te = new TwistEdge( **strand, vtx, strandRadius, m_steppers.back() );

            if( vtx > 0 ){
                te->prev = last;
                te->prev->next = te;
            }
            originalProxies.push_back( te );
            last = te;
        }
    }

    for ( auto controller = m_meshScriptingControllers.begin(); controller != m_meshScriptingControllers.end(); ++controller )
    {
        const TriangularMesh* const mesh = ( *controller )->getCurrentMesh();

        for ( unsigned f = 0; f < mesh->nf(); ++f ){
            originalProxies.push_back( new FaceProxy( f, *controller ) );
        }
    }

}


bool StrandImplicitManager::isCollisionInvariantCT( const Scalar dt )
{
    unsigned prevNumCollisions = m_mutualContacts.size();
    m_collisionDetector->clear();
    m_collisionDetector->buildBVH( false );

    EdgeFaceIntersection::s_doProximityDetection = true;
    m_collisionDetector->findCollisions( !m_params.m_useCTRodRodCollisions , false , true ); // set whether cd should do ctc for hairs, etc
    doContinuousTimeDetection( dt );

    bool collisionFreeCT = true;

    if( m_mutualContacts.size() == prevNumCollisions ){ 
        return true;
    }
    else{
        std::list<CollisionBase*>& collisionsList = m_collisionDetector->getContinuousTimeCollisions();
        for ( auto collIt = collisionsList.begin(); collIt != collisionsList.end(); ++collIt )
        {
            EdgeEdgeCollision* const ctCollision = dynamic_cast<EdgeEdgeCollision*>( *collIt );
            if ( !ctCollision ){
                continue;
            }
            if( ctCollision->m_firstStrand->getGlobalIndex() != ctCollision->m_secondStrand->getGlobalIndex() )
            {
                // ctCollision->printShort( std::cout ); cout << endl;
                collisionFreeCT = false;
                break;
            }
        }
    }

    return collisionFreeCT;
}



void StrandImplicitManager::setupContinuousTimeCollisions()
{
    unsigned numTunneledBands;
    do{
        m_collisionDetector->m_proxyHistory->repeatCD = false;

        m_collisionDetector->buildBVH( false );
        EdgeFaceIntersection::s_doProximityDetection = true;
        m_collisionDetector->findCollisions( !m_params.m_useCTRodRodCollisions , false , false ); // set whether cd should do ctc for hairs, etc

    } while( m_collisionDetector->m_proxyHistory->repeatCD );
}

void StrandImplicitManager::doProximityMeshHairDetection( Scalar dt )
{
    std::list<CollisionBase*>& collisionsList = m_collisionDetector->getProximityCollisions();

    unsigned nInt = 0;
    for ( auto intersection = collisionsList.begin(); intersection != collisionsList.end(); ++intersection )
    {
        if ( EdgeFaceIntersection* efi = dynamic_cast<EdgeFaceIntersection*>( *intersection ) )
        {
            EdgeProxy* edge = efi->getEdge();

            const unsigned strIdx = edge->getStrand().getGlobalIndex();
            const unsigned edgeIdx = edge->getVertexIndex();

            ProximityCollision edgeFaceCollision;

            edgeFaceCollision.normal = efi->getFaceNormal();
            edgeFaceCollision.mu = sqrt( efi->faceFrictionCoefficient() * edge->getStrand().collisionParameters().frictionCoefficient( edgeIdx ) );
            edgeFaceCollision.adhesion = ( efi->faceAdhesionCoefficient() + edge->getStrand().collisionParameters().adhesionCoefficient( edgeIdx ) ) / dt;

            edgeFaceCollision.objects.second.globalIndex = -1;
            edgeFaceCollision.objects.second.vertex = efi->getFaceId();
            edgeFaceCollision.objects.second.freeVel = efi->getFaceVelocity( dt ); //+ ContinuousTimeCollision::s_extraRadius * collision.normal ;

            if ( addExternalContact( strIdx, edgeIdx, efi->getEdgeAbscissa(), edgeFaceCollision ) ){
                ++nInt;
            }
        }
    }
}


void StrandImplicitManager::doContinuousTimeDetection( Scalar dt )
{
    m_num_ct_hair_hair_col = 0;
    std::list<CollisionBase*>& collisionsList = m_collisionDetector->getContinuousTimeCollisions();
    if ( collisionsList.empty() && m_collisionDetector->m_proxyHistory->tunnelingBands.size() == 0 ){
        return;
    }

    // In attempt to eliminate duplicates
    collisionsList.sort( compareCT );

    unsigned nInt = 0;
    unsigned nCTCD = 0;

    CollisionBase *previous = NULL;
    for ( auto collIt = collisionsList.begin(); collIt != collisionsList.end(); ++collIt )
    {
        if ( previous && !compareCT( previous, *collIt ) ) // Therefore they are equal
        {
            continue;
        }
        previous = *collIt;

        ContinuousTimeCollision* const ctCollision = dynamic_cast<ContinuousTimeCollision*>( *collIt );
        if ( !ctCollision ) continue;

        ProximityCollision collision;
        ElasticStrand* const strand = ctCollision->getFirstStrand();
        const unsigned edgeIdx = ctCollision->getFirstVertex();

        collision.m_originalCTCollision = ctCollision;
        collision.normal = ctCollision->normal();

        if ( FaceCollision* fCollision = dynamic_cast<FaceCollision*>( ctCollision ) ){
            collision.mu = sqrt( fCollision->faceFrictionCoefficient() * strand->collisionParameters().frictionCoefficient( edgeIdx ) );
            collision.adhesion = ( fCollision->faceAdhesionCoefficient() + strand->collisionParameters().adhesionCoefficient( edgeIdx ) ) / dt;
        }
        else if ( m_params.m_useCTRodRodCollisions )
        {
            EdgeEdgeCollision* eCollision = dynamic_cast<EdgeEdgeCollision*>( ctCollision );
            if ( !eCollision ) continue;
                    
            ElasticStrand* const strand1 = eCollision->getFirstStrand() ;
            const unsigned edgeIdx1 = eCollision->getFirstVertex() ;
            ElasticStrand* const strand2 = eCollision->getSecondStrand() ;
            const unsigned edgeIdx2 = eCollision->getSecondVertex() ;
            
            collision.distance = 0. ;
            collision.mu = sqrt( strand1->collisionParameters().frictionCoefficient( edgeIdx1 ) * strand2->collisionParameters().frictionCoefficient( edgeIdx2 ) ) ;
            collision.adhesion = ( strand1->collisionParameters().adhesionCoefficient( edgeIdx1 ) + strand2->collisionParameters().adhesionCoefficient( edgeIdx2 ) ) / dt ;
            
            collision.objects.first.globalIndex = strand1->getGlobalIndex() ;
            collision.objects.first.vertex = edgeIdx1 ;

            // Scalar multp = 1.0 + hIter * 0.1;

            // collision.objects.first.freeVel = - ( eCollision->getPenetrationVector() );
            // collision.objects.first.freeVel *= multp;

            // cout << "first freeVel reg: " << collision.objects.first.freeVel << endl;
            collision.objects.second.globalIndex = strand2->getGlobalIndex() ;
            collision.objects.second.vertex = edgeIdx2 ;
            // collision.objects.second.freeVel =  ( eCollision->getPenetrationVector() );
            // collision.objects.second.freeVel *= multp;


            collision.objects.first.freeVel = Vec3x::Zero();
            collision.objects.second.freeVel = Vec3x::Zero();            

            collision.objects.first.abscissa = eCollision->getFirstAbscissa();
            collision.objects.second.abscissa = eCollision->getSecondAbscissa();

            ++m_num_ct_hair_hair_col;
            unsigned outof = std::floor( 100. / m_params.m_percentCTRodRodCollisionsAccept );
#pragma omp critical
{
            if( !( m_num_ct_hair_hair_col % outof ) ) // for sample testing
            {
                collision.swapIfNecessary();

                int redundant = redundantConstraintIndex( collision, currentStepCollisions );
                if( redundant != -1 )
                {
                    collision.objects.first.freeVel *= 1.0;
                    collision.objects.second.freeVel *= 1.0;
                }
                else
                {
                    currentStepCollisions.push_back( collision );
                }

                m_mutualContacts.push_back( collision );
                Vec3x pos = strand1->getVertex( edgeIdx1 );
                m_steppers[0]->m_strand.rod_data->m_renderer->addVert( pos );
                pos = strand2->getVertex( edgeIdx2 );
                m_steppers[0]->m_strand.rod_data->m_renderer->addVert( pos );
                // cout << "CTCD [" << m_mutualContacts.size() << "] "; collision.printShort( std::cout ); cout << endl;
                ++nCTCD;
            }
}            
            continue;
        }


        collision.objects.second.globalIndex = -1;
        const unsigned strIdx = strand->getGlobalIndex();
        const Vec3x offset = ctCollision->offset() / dt;

        if ( VertexFaceCollision* vfCollision = dynamic_cast<VertexFaceCollision*>( ctCollision ) )
        {
            collision.objects.second.vertex = vfCollision->face()->uniqueId();
            collision.objects.second.freeVel = vfCollision->meshVelocity( dt ) + offset;
            if ( addExternalContact( strIdx, edgeIdx, 0, collision ) ){
                ++nInt;
            }
        }
        else if ( EdgeFaceCollision* efCollision = dynamic_cast<EdgeFaceCollision*>( ctCollision ) )
        {
            collision.objects.second.vertex = efCollision->faceEdgeId();
            collision.objects.second.freeVel = efCollision->meshVelocity( dt ) + offset;
            if ( addExternalContact( strIdx, edgeIdx, efCollision->abscissa(), collision ) ){
                ++nInt;
            }
        }

    }


    std::cout << "CTCD " << nCTCD << std::endl;

}

bool StrandImplicitManager::addExternalContact( const unsigned strIdx, const unsigned edgeIdx, const Scalar abscissa, const ProximityCollision& externalContact )
{
    if ( !acceptsCollision( *m_strands[strIdx], edgeIdx, abscissa, externalContact.normal ) )
        return false;

    // Discard collisions if their normal is almost parallel to the strand edge and would cause stretching
    const int prevEdge = abscissa == 0. ? edgeIdx - 1 : edgeIdx;
    const Vec3x edge = m_strands[strIdx]->getCurrentTangent( prevEdge );
    if ( edge.dot( externalContact.normal ) > ALMOST_PARALLEL_COS ){
        std::cout << "CULLING out external contact near at edgeIdx: " << edgeIdx << " because ALMOST PARALLEL COS " << std::endl;
        return false;
    }

    m_externalContacts[strIdx].push_back( externalContact );
    ProximityCollision &collision = m_externalContacts[strIdx].back();

    collision.objects.first.globalIndex = strIdx;
    collision.objects.first.vertex = edgeIdx;
    collision.objects.first.abscissa = abscissa;
    collision.objects.first.freeVel.setZero();

    return true;
}

bool StrandImplicitManager::acceptsCollision( const ElasticStrand& strand, int edgeIdx, Scalar localAbscissa, const Vec3x & )
{
    if ( !edgeIdx || ( edgeIdx == 1 && localAbscissa < SECOND_EDGE_MIN_CONTACT_ABSCISSA ) ){
        return false;
    }
    if ( strand.getCurvilinearAbscissa( edgeIdx, localAbscissa ) < strand.collisionParameters().rootImmunityLength() ){
        return false;
    }
    return true;
}

void StrandImplicitManager::setupProximityCollisions( Scalar dt )
{ // [H] this is purely Proximity CD, using spatial hash

    // do hair-hair setup basics and proceed to proximity collisions if requested
    if ( !m_params.m_useProxRodRodCollisions ) return;

    if ( !m_hashMap ){
        m_hashMap = new SpatialHashMapT();
    }

    // spatial grid size determined by max of avg rod radius and max edge length (is this what we want?)
    // Compute maximum number of vertices and min colliding radius
    unsigned maxNumVert = 0;
    Scalar meanRadius = 0.;
    Scalar maxEdgeLen = 0.;
    for ( unsigned i = 0; i < m_strands.size(); ++i )
    {
        const unsigned nv = m_strands[i]->getNumVertices();
        if ( nv > maxNumVert ){
            maxNumVert = nv;
        }

        const Scalar rootRad = m_strands[i]->collisionParameters().selfCollisionsRadius( 0 );
        meanRadius += rootRad;
        maxEdgeLen = std::max( maxEdgeLen, m_strands[i]->getTotalRestLength() / nv );
    }
    m_hashMap->setCellSize( std::max( maxEdgeLen, meanRadius / m_strands.size() ) );

    m_hashMap->batchUpdate( m_strands, maxNumVert );
        
    // compute will actually do nothing, except initializing result
    // The collisions will actually be retrieved by batches of 10 objects at each call to result.next ([H], what?)
    SpatialHashMapT::Result result( true, 10 );
    m_hashMap->compute( result );

    unsigned nRough = 0, nExt = 0;
    unsigned nProx = 0;

    // Processes batches of collision in parallel
    SpatialHashMapT::Result::Collisions collisions;
    
#pragma omp parallel private( collisions ) reduction( + : nRough, nExt )
    while ( result.next( collisions ) )
    {

        const unsigned nColObj = collisions.size();

        std::vector<SpatialHashMapT::Result::Collisions::const_iterator> iters( nColObj );

        unsigned i = 0;
        for ( SpatialHashMapT::Result::Collisions::const_iterator first = collisions.begin(); first != collisions.end(); ++first )
        {
            iters[i++] = first;
        }

        for ( unsigned itIdx = 0; itIdx < nColObj; ++itIdx )
        {
            SpatialHashMapT::Result::Collisions::const_iterator &first = iters[itIdx];
            ElasticStrand* sP = first->first;
            const CollisionParameters &cpP = sP->collisionParameters();

            if ( !cpP.createsSelfCollisions() )
                continue;

            for ( auto second = first->second.begin(); second != first->second.end(); ++second )
            {
                ElasticStrand* sQ = second->first;

                const CollisionParameters &cpQ = sQ->collisionParameters();
                if ( !cpQ.createsSelfCollisions() )
                    continue;

                for ( auto collision = second->second.begin(); collision != second->second.end();
                        ++collision )
                {
                    int iP = collision->first;
                    int iQ = collision->second;

                    if ( ( sP == sQ && std::abs( iP - iQ ) < 4 ) || !( iP && iQ ) )
                        continue;

                    ++nRough;

                    // [H] Narrow detection phase:
                    Vec3x normal;
                    Scalar s, t, d;
                    if ( !analyseRoughRodRodCollision( sP, sQ, iP, iQ, normal, s, t, d ) ){
                        continue;
                    }

                    bool acceptFirst  = cpP.reactsToSelfCollisions() && acceptsCollision( *sP, iP, s, normal );
                    bool acceptSecond = cpQ.reactsToSelfCollisions() && acceptsCollision( *sQ, iQ, t, -normal );

                    // && -> Discard collisions on root immunity length
                    // || -> make them as external objects
                    if ( !( acceptFirst || acceptSecond ) )
                        continue;
                    if ( cpP.usesFakeLayering() && cpQ.usesFakeLayering() )
                    {
                        if ( !( acceptFirst && acceptSecond ) )
                            continue;
                    }

                    ProximityCollision mutualContact;
                    mutualContact.normal = normal;
                    mutualContact.distance = d;
                    mutualContact.mu = sqrt( cpP.frictionCoefficient( iP ) * cpQ.frictionCoefficient( iQ ) );
                    mutualContact.adhesion = ( cpP.adhesionCoefficient( iP ) + cpQ.adhesionCoefficient( iQ ) ) / dt;

                    mutualContact.objects.first.globalIndex = sP->getGlobalIndex();
                    mutualContact.objects.first.vertex = iP;
                    mutualContact.objects.first.abscissa = s;
                    mutualContact.objects.first.freeVel.setZero();

                    mutualContact.objects.second.globalIndex = sQ->getGlobalIndex();
                    mutualContact.objects.second.vertex = iQ;
                    mutualContact.objects.second.abscissa = t;
                    mutualContact.objects.second.freeVel.setZero();

                    if ( acceptFirst && acceptSecond )
                    {
                        mutualContact.swapIfNecessary();
#pragma omp critical
                        {
                            // cout << "PROX "; mutualContact.printShort( std::cout ); cout << endl;
                            m_mutualContacts.push_back( mutualContact );
                            Vec3x pos = sP->getVertex( iP );
                            m_steppers[0]->m_strand.rod_data->m_renderer->addVert( pos );                            
                            ++nProx;
                        }
                    }
                    else
                    {
                        ++nExt;
                        makeExternalContact( mutualContact, acceptFirst );
                    }
                }
            }
        }
    }

    if ( m_params.m_simulationManager_limitedMemory )
    {
        delete m_hashMap;
        m_hashMap = NULL;
    }
    else
    {
        m_hashMap->clear();
    }

    std::cout << "Prox " << nProx << std::endl;
}

void StrandImplicitManager::makeExternalContact( ProximityCollision& externalContact, bool onFirstObject )
{
    int extObjId;
    if ( onFirstObject )
    {
        extObjId = externalContact.objects.second.globalIndex;
        externalContact.objects.second.globalIndex = -1;
    }
    else
    {
        extObjId = externalContact.objects.first.globalIndex;
        externalContact.objects.first.globalIndex = -1;
    }

    // This will put the external object on c.objects.second
    externalContact.swapIfNecessary();

    const VecXx &velocities = m_steppers[extObjId]->velocities();
    const int iQ = externalContact.objects.second.vertex;
    externalContact.objects.second.vertex = -extObjId;
    const Scalar t = externalContact.objects.second.abscissa;

    externalContact.objects.second.freeVel = ( 1 - t ) * velocities.segment<3>( 4 * iQ ) + t * velocities.segment<3>( 4 * iQ + 4 );

#pragma omp critical
    {
        m_externalContacts[externalContact.objects.first.globalIndex].push_back( externalContact );
    }
}

void StrandImplicitManager::computeDeformationGradient( ProximityCollision::Object &object ) const
{
    m_strands[object.globalIndex]->getFutureState().computeDeformationGradient( 
        object.vertex, object.abscissa, m_steppers[object.globalIndex]->velocities(), object.defGrad, object.freeVel );
}

void StrandImplicitManager::setupDeformationBasis( ProximityCollision &collision ) const
{
    const ProximityCollision* oldCollision = m_collisionDatabase.find( collision );
    if ( oldCollision )
    {
        collision.force = oldCollision->force; // Warm start
        collision.updateTransformationMatrix( oldCollision->transformationMatrix );
    }
    else{
        collision.generateTransformationMatrix();
    }

    computeDeformationGradient( collision.objects.first );
    if ( collision.objects.second.globalIndex != -1 ){
        computeDeformationGradient( collision.objects.second );
    }
}




bool StrandImplicitManager::needsExternalSolve( unsigned strandIdx ) const
{
    return !m_externalContacts[strandIdx].empty();
}


void StrandImplicitManager::solveSingleObject( unsigned objectIdx, bool asFailSafe, bool nonLinear )
{
    CollidingGroup collisionGroup;
    collisionGroup.first[objectIdx] = 0;

    for ( unsigned i = 0; i < m_mutualContacts.size(); ++i )
    {
        ProximityCollision &mutualCollision = m_mutualContacts[i];

        setupDeformationBasis( mutualCollision );

        const unsigned s1 = mutualCollision.objects.first.globalIndex;
        // const unsigned s2 = mutualCollision.objects.second.globalIndex;



        // if( s1 == objectIdx )
        // {
        //     collisionGroup.second.push_back( mutualCollision );
        // }
    }

    solveCollidingGroup( collisionGroup, asFailSafe, nonLinear );
}

// DK: another place the failsafes kick in -- here we have it outside the GS solve *LOOK HERE*
void StrandImplicitManager::solveCollidingGroup( CollidingGroup& collisionGroup, bool asFailSafe, bool nonLinear )
{
    if ( collisionGroup.first.empty() ) return;

    std::vector<unsigned> globalIds;
    std::vector<ProximityCollision*> colPointers;

    bool accept = false;
    {
        if( useQL ){
            accept = solveQLFrictionProblem( collisionGroup, globalIds, colPointers );
        }
        else{
            VecXx vels;
            VecXx impulses;
            bogus::MecheFrictionProblem mecheProblem;
            VecXu startDofs;
            VecXu nDofs;
            int numSubSystems;
            if( assembleBogusFrictionProblem( collisionGroup, mecheProblem, globalIds, colPointers, vels, impulses, startDofs, nDofs, numSubSystems ) )
            {
                // if( asFailSafe ){
                // cout << "bogus attempted" << endl;
                    accept = solveBogusFrictionProblem( mecheProblem, globalIds, asFailSafe, nonLinear, vels, impulses, numSubSystems );
                postProcessBogusFrictionProblem( accept, collisionGroup, mecheProblem, globalIds, colPointers, vels, impulses, startDofs, nDofs);
                // cout << "accepted: " << accept << endl;           
                // }

            }
        }
    }

    /* If either the solver result was really bad or one strand was stretching,
     we trigger a hierarchy of failsafes.

     - First retry with the non-linear solver if this was not already the case
     - Then discard all mutual contacts, solve each strand with non-linear solver
     - Finally solve each strand with linear solver but with length constraints
     - If that is still failing, trust the ImplicitStepper's reprojection
     */

    // DK: failsafes here:
    bool mustRetry = !accept;
    for ( unsigned i = 0; i < globalIds.size(); ++i )
    {
        const unsigned sIdx = globalIds[i];
        if ( m_steppers[sIdx]->lastStepWasRejected() )
        {
            cout << "lastStepWasRejected... " << endl;
            mustRetry = true;
            break;
        }
    }
/*
    if ( mustRetry )
    {
        std::cout << "Collision routine failed, must retry" << std::endl;

        if ( globalIds.size() > 1 )
        {
            // Failed again, drop mutual collisions
#pragma omp parallel for
            for ( unsigned i = 0; i < globalIds.size(); ++i )
            {
                const unsigned sIdx = globalIds[i];
                if ( !accept || ( m_steppers[sIdx]->lastStepWasRejected() ) ){
                    cout << "solveSingleObject called with true" << endl;
                    solveSingleObject( globalIds[i], true, m_params.m_useNonLinearAsFailsafe || m_params.m_alwaysUseNonLinear );
                }

            }
        }
        // DK: (length projection) otherwise if nonlinear "always" try length-constraint as failsafe...
        else if ( m_params.m_useLengthConstraints && ( !asFailSafe || nonLinear ) )
        {
            ImplicitStepper& stepper = *m_steppers[globalIds[0]];
            std::cerr << "StrandImplicitManager::solveCollidingGroup, enabling BilateralConstraints... not yet supported" << std::endl;

            // Single object failed, retry with length constraints
            stepper.enableBilateralConstraints();
            solveCollidingGroup( collisionGroup, true, false );
        }

    }
    */

    if ( !asFailSafe )
    {

#pragma omp parallel for
        for ( unsigned i = 0; i < colPointers.size(); ++i )
        {
            ProximityCollision& col = *colPointers[i];
            delete col.objects.first.defGrad;
            col.objects.first.defGrad = NULL;
            if ( col.objects.second.globalIndex != -1 )
            {
                delete col.objects.second.defGrad;
                col.objects.second.defGrad = NULL;
            }
        }

    }
}


struct ProxColPointerCompare
{
    bool operator()( const ProximityCollision* lhs, const ProximityCollision* rhs ) const
    {
        return *lhs < *rhs;
    }
};

void StrandImplicitManager::computeCollidingGroups( const ProximityCollisions &origMutualCollisions, Scalar  )
{
    if( origMutualCollisions.empty() ){
        // cout <<"[StrandImplicitManager::computeCollidingGroups] no collisions detected " <<endl;
        return;
    }

    ProximityCollisions mutualCollisions;
    mutualCollisions.reserve( origMutualCollisions.size() );

    if ( m_params.m_pruneSelfCollisions )
    {
        pruneCollisions( origMutualCollisions, mutualCollisions, m_params.m_stochasticPruning );
    }
    else
    {
        for ( unsigned i = 0; i < origMutualCollisions.size(); ++i )
        {
            const ProximityCollision &proxyCol = origMutualCollisions[i];
            const unsigned s1 = proxyCol.objects.first.globalIndex;
            const unsigned s2 = proxyCol.objects.second.globalIndex;

            if ( m_steppers[s1]->refusesMutualContacts()
                    || m_steppers[s2]->refusesMutualContacts() )
            {
                ProximityCollision copy( proxyCol );
                makeExternalContact( copy, m_steppers[s2]->refusesMutualContacts() );
            }
            else
            {
                mutualCollisions.push_back( proxyCol );
            }
        }
        if ( m_params.m_useDeterministicSolver )
        {
            std::sort( mutualCollisions.begin(), mutualCollisions.end() );
        }
    }

    m_statMutualCollisions += mutualCollisions.size();

    // For each strand, list all other contacting ones
    std::vector<std::deque<unsigned> > objsGroups( m_strands.size() );
    for ( unsigned i = 0; i < mutualCollisions.size(); ++i )
    {
        const ProximityCollision &proxyCol = mutualCollisions[i];
        const unsigned s1 = proxyCol.objects.first.globalIndex;
        const unsigned s2 = proxyCol.objects.second.globalIndex;

        objsGroups[s1].push_back( s2 );
        objsGroups[s2].push_back( s1 );
    }

    // Extract connected subgraphs from the global constraint graph
    for ( unsigned s1 = 0; s1 < objsGroups.size(); ++s1 )
    {
        if ( m_collidingGroupsIdx[s1] != -1 || objsGroups[s1].empty() )
            continue;

        const unsigned groupIdx = m_collidingGroups.size();

        m_collidingGroups.push_back( CollidingGroup() );
        CollidingGroup &cg = m_collidingGroups.back();

        m_collidingGroupsIdx[s1] = groupIdx;
        cg.first[s1] = 0;

        std::deque<unsigned> toVisit = objsGroups[s1];

        while ( toVisit.size() )
        {
            const unsigned s2 = toVisit.front();
            toVisit.pop_front();

            if ( m_collidingGroupsIdx[s2] != -1 )
                continue;

            m_collidingGroupsIdx[s2] = groupIdx;
            cg.first[s2] = 0;

            toVisit.insert( toVisit.end(), objsGroups[s2].begin(), objsGroups[s2].end() );
        }

    }

    for ( unsigned i = 0; i < mutualCollisions.size(); ++i )
    {
        const ProximityCollision &mutualCollision = mutualCollisions[i];
        const unsigned s1 = mutualCollision.objects.first.globalIndex;

        m_collidingGroups[m_collidingGroupsIdx[s1]].second.push_back( mutualCollision );
    }

    // Index in group needs to be the same that when using std::map iterator
#pragma omp parallel for
    for ( unsigned i = 0; i < m_collidingGroups.size(); ++i )
    {
        unsigned k = 0;
        IndicesMap &indices = m_collidingGroups[i].first;
        for ( IndicesMap::iterator it = indices.begin(); it != indices.end(); ++it )
        {
            it->second = k++;
        }
    }
}
