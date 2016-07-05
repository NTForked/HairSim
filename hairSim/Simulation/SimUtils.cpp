#include "Simulation.h"
#include "../Collision/Collision.h"
#include "../Collision/ElementProxy.h"
#include "../Collision/CollisionDetector.h"
#include "../Collision/CollisionUtils/SpatialHashMap.hh"
#include "../Collision/CollisionUtils/CollisionUtils.h"
#include "../Collision/VertexFaceCollision.h"
#include "../Collision/EdgeFaceCollision.h"
#include <Eigen/Sparse>

#define SECOND_EDGE_MIN_CONTACT_ABSCISSA 0.0001
#define ALMOST_PARALLEL_COS 0.96592582628 // cos( Pi/12 )

void Simulation::accumulateProxies( std::vector< ElementProxy* >& originalProxies, const std::vector< TriMesh* >& meshes )
{
    // Create the stepper for each strand
    // Also store edge proxies for collision detection
    
    TwistEdge* last = NULL;
    for( std::vector<ElasticStrand*>::const_iterator strand = m_strands.begin(); strand != m_strands.end(); ++strand )
    {
        m_steppers.push_back( new ImplicitStepper( **strand, m_params ) );

        for ( int vtx = 0; vtx < ( *strand )->getNumEdges(); ++vtx )
        {
            TwistEdge* te = new TwistEdge( **strand, vtx );

            if( vtx > 0 ){
                te->prev = last;
                te->prev->next = te;
            }
            originalProxies.push_back( te );
            last = te;
        }
    }

    for( unsigned m = 0; m < meshes.size(); ++m )
    {
        TriMesh* mesh = meshes[m];
        for( unsigned f = 0; f < mesh->nf(); ++f ){
            originalProxies.push_back( new FaceProxy( f, mesh->controller() ) );
        }
    }
}

bool Simulation::isCollisionInvariantCT( const Scalar dt )
{
    bool collisionFreeCT = true;

    unsigned prevNumCollisions = m_mutualContacts.size();
    { // Detect and Preprocess Collisions
        m_collisionDetector->clear();
        m_collisionDetector->buildBVH( false );
        m_collisionDetector->findCollisions( false );
        preProcessContinuousTimeCollisions( dt );
    }

    std::cerr << "this may have to change to just looking for number greater than 0" << std::endl;
    if( m_mutualContacts.size() == prevNumCollisions ){ 
        return true;
    }
    else{
        const std::list< Collision* >& collisionsList = m_collisionDetector->getCollisions();
        for ( auto collIt = collisionsList.begin(); collIt != collisionsList.end(); ++collIt )
        {
            EdgeEdgeCollision* const ctCollision = dynamic_cast<EdgeEdgeCollision*>( *collIt );
            if ( !ctCollision ){
                continue;
            }
            collisionFreeCT = false;
            break;
        }
    }

    return collisionFreeCT;
}

void Simulation::gatherProximityRodRodCollisions( Scalar dt )
{ // Detect with SpatialHash, then preprocess

    if( !m_params.m_useProxRodRodCollisions ) return;

    if( !m_hashMap ){
        m_hashMap = new SpatialHashMapT();
    }

    // spatial grid size determined by max of ( avg rod radius and max edge length )
    // Compute maximum number of vertices and min colliding radius
    unsigned maxNumVert = 0;
    Scalar meanRadius = 0.;
    Scalar maxEdgeLen = 0.;
    for ( unsigned i = 0; i < m_strands.size(); ++i )
    {
        const unsigned nv = m_strands[i]->getNumVertices();
        if( nv > maxNumVert ){
            maxNumVert = nv;
        }

        const Scalar rootRad = m_strands[i]->collisionParameters().getLargerCollisionsRadius( 0 );
        meanRadius += rootRad;
        maxEdgeLen = std::max( maxEdgeLen, m_strands[i]->getTotalRestLength() / nv );
    }
    m_hashMap->setCellSize( std::max( maxEdgeLen, meanRadius / m_strands.size() ) );
    m_hashMap->batchUpdate( m_strands, maxNumVert );
        
    // compute will actually do nothing, except initializing result
    // The collisions will actually be retrieved by batches of 10 objects at each call to result.next
    SpatialHashMapT::Result result( true, 10 );
    m_hashMap->compute( result );

    unsigned nRough = 0, nExt = 0, nProx = 0;

    // Processes batches of collision in parallel
    SpatialHashMapT::Result::Collisions collisions;
    
#pragma omp parallel private( collisions ) reduction( + : nRough, nExt )
    while( result.next( collisions ) )
    {
        const unsigned nColObj = collisions.size();
        std::vector<SpatialHashMapT::Result::Collisions::const_iterator> iters( nColObj );

        unsigned i = 0;
        for( SpatialHashMapT::Result::Collisions::const_iterator first = collisions.begin(); first != collisions.end(); ++first )
        {
            iters[i++] = first;
        }

        for( unsigned itIdx = 0; itIdx < nColObj; ++itIdx )
        {
            SpatialHashMapT::Result::Collisions::const_iterator &first = iters[itIdx];
            ElasticStrand* sP = first->first;
            const CollisionParameters &cpP = sP->collisionParameters();

            for( auto second = first->second.begin(); second != first->second.end(); ++second )
            {
                ElasticStrand* sQ = second->first;
                const CollisionParameters &cpQ = sQ->collisionParameters();

                for( auto collision = second->second.begin(); collision != second->second.end(); ++collision )
                {
                    int iP = collision->first;
                    int iQ = collision->second;

                    if( ( sP == sQ && std::abs( iP - iQ ) < 4 ) || !( iP && iQ ) ){
                        continue;
                    }

                    ++nRough;

                    // [H] Narrow detection phase:
                    Vec3 normal;
                    Scalar s, t, d;
                    if( !analyseRoughRodRodCollision( sP, sQ, iP, iQ, normal, s, t, d ) ){
                        continue;
                    }

                    bool acceptFirst  = acceptsCollision( *sP, iP, s );
                    bool acceptSecond = acceptsCollision( *sQ, iQ, t );

                    if( !acceptFirst && !acceptSecond ){
                        continue;
                    }

                    CollidingPair mutualContact;
                    mutualContact.m_normal = normal;
                    mutualContact.m_mu = sqrt( cpP.m_frictionCoefficient * cpQ.m_frictionCoefficient );

                    mutualContact.objects.first.globalIndex = sP->getGlobalIndex();
                    mutualContact.objects.first.vertex = iP;
                    mutualContact.objects.first.abscissa = s;
                    mutualContact.objects.first.worldVel.setZero();

                    mutualContact.objects.second.globalIndex = sQ->getGlobalIndex();
                    mutualContact.objects.second.vertex = iQ;
                    mutualContact.objects.second.abscissa = t;
                    mutualContact.objects.second.worldVel.setZero();

                    if( acceptFirst && acceptSecond )
                    {
                        mutualContact.swapIfNecessary();
#pragma omp critical
                        {
                            m_mutualContacts.push_back( mutualContact );
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
    std::cout << "Prox " << nProx << std::endl;

    if( m_params.m_simulationManager_limitedMemory )
    {
        delete m_hashMap;
        m_hashMap = NULL;
    }
    else
    {
        m_hashMap->clear();
    }
}

bool Simulation::acceptsCollision( const ElasticStrand& strand, int edgeIdx, Scalar localAbscissa )
{
    if( !edgeIdx || ( edgeIdx == 1 && localAbscissa < SECOND_EDGE_MIN_CONTACT_ABSCISSA ) )
    { // very close to root, ignore
        return false;
    }
    if( strand.getCurvilinearAbscissa( edgeIdx, localAbscissa ) < strand.collisionParameters().m_rootImmunityLength ){
        return false;
    }
    return true;
}

void Simulation::makeExternalContact( CollidingPair& externalContact, bool onFirstObject )
{
    int extObjId;
    if( onFirstObject )
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

    const VecXx &velocities = m_steppers[extObjId]->m_futureVelocities;
    const int iQ = externalContact.objects.second.vertex;
    externalContact.objects.second.vertex = -extObjId;
    const Scalar t = externalContact.objects.second.abscissa;

    // free vel is just the velocity of the second object, not relative velocity
    externalContact.objects.second.worldVel = ( 1 - t ) * velocities.segment<3>( 4 * iQ ) + t * velocities.segment<3>( 4 * iQ + 4 );

#pragma omp critical
    {
        m_externalContacts[ externalContact.objects.first.globalIndex ].push_back( externalContact );
    }
}

void Simulation::detectContinuousTimeCollisions()
{ // dont clear collisions here since we only process them after gathering all of them in this loop
    do{ 
        m_collisionDetector->m_proxyHistory->repeatCD = false;
        m_collisionDetector->buildBVH( false );
        m_collisionDetector->findCollisions( false );
    } while( m_collisionDetector->m_proxyHistory->repeatCD );
}

void Simulation::preProcessContinuousTimeCollisions( Scalar dt )
{
    std::list< Collision* >& collisionsList = m_collisionDetector->getCollisions();
    if( collisionsList.empty() ){
        return;
    }
    collisionsList.sort( compareCT );

    unsigned nInt = 0; // External
    unsigned nCTCD = 0; // RodRod
    Collision *previous = NULL;
    for( auto collIt = collisionsList.begin(); collIt != collisionsList.end(); ++collIt )
    {
        if( previous && !compareCT( previous, *collIt ) ) // Therefore they are equal
        {
            continue;
        }
        previous = *collIt;

        Collision* const ctCollision = dynamic_cast< Collision* >( *collIt );
        if( !ctCollision ) continue;

        CollidingPair collision;
        ElasticStrand* const strand1 = ctCollision->getFirstEdgeProxy()->getStrandPointer();
        const unsigned edgeIdx1 = ctCollision->getFirstEdgeProxy()->getVertexIndex();

        collision.m_normal = ctCollision->normal();

        collision.objects.first.globalIndex = strand1->getGlobalIndex();
        collision.objects.first.vertex = edgeIdx1;
        EdgeEdgeCollision* eCollision = dynamic_cast<EdgeEdgeCollision*>( ctCollision );

        if( m_params.m_useCTRodRodCollisions && eCollision != NULL  )
        {

            collision.objects.first.abscissa = eCollision->getFirstAbscissa();

            collision.objects.first.worldVel = 
                (1.0 - eCollision->getFirstAbscissa() ) * m_steppers[strand1->getGlobalIndex()]->m_futureVelocities.segment<3>( 4 * edgeIdx1 ) 
                + eCollision->getFirstAbscissa() * m_steppers[strand1->getGlobalIndex()]->m_futureVelocities.segment<3>( 4 * (edgeIdx1 + 1) ); //not relative vel, but just vel of each

            ElasticStrand* const strand2 = eCollision->getSecondEdgeProxy()->getStrandPointer();
            const unsigned edgeIdx2 = eCollision->getSecondEdgeProxy()->getVertexIndex();
            
            collision.m_mu = sqrt( strand1->collisionParameters().m_frictionCoefficient * strand2->collisionParameters().m_frictionCoefficient ) ;
            
            collision.objects.second.globalIndex = strand2->getGlobalIndex();
            collision.objects.second.vertex = edgeIdx2;
            collision.objects.second.abscissa = eCollision->getSecondAbscissa();
            collision.objects.second.worldVel = 
                (1.0 - eCollision->getSecondAbscissa() ) * m_steppers[strand2->getGlobalIndex()]->m_futureVelocities.segment<3>( 4 * edgeIdx2 ) 
                + eCollision->getSecondAbscissa() * m_steppers[strand2->getGlobalIndex()]->m_futureVelocities.segment<3>( 4 * (edgeIdx2 + 1) ); //not relative vel, but just vel of each

            collision.swapIfNecessary();

#pragma omp critical
            {
                m_mutualContacts.push_back( collision );
                ++nCTCD;
            }            
            continue;
        }

        collision.objects.second.globalIndex = -1;
        if( VertexFaceCollision* vfCollision = dynamic_cast<VertexFaceCollision*>( ctCollision ) )
        {
            collision.objects.second.abscissa = 0;

            const Vec3 offset = vfCollision->offset( ctCollision->normal() ) / dt;            
            collision.m_mu = sqrt( vfCollision->faceFrictionCoefficient() * strand1->collisionParameters().m_frictionCoefficient );
            collision.objects.second.vertex = vfCollision->face()->uniqueId();
            collision.objects.second.worldVel = vfCollision->meshVelocity( dt ) + offset;
            if( addExternalContact( strand1->getGlobalIndex(), edgeIdx1, 0, collision ) ){
                ++nInt;
            }
        }
        else if( EdgeFaceCollision* efCollision = dynamic_cast<EdgeFaceCollision*>( ctCollision ) )
        {
            collision.objects.first.abscissa = efCollision->abscissa();

            const Vec3 offset = efCollision->offset( ctCollision->normal() ) / dt;
            collision.m_mu = sqrt( efCollision->faceFrictionCoefficient() * strand1->collisionParameters().m_frictionCoefficient );
            collision.objects.second.vertex = efCollision->faceEdgeId();
            collision.objects.second.worldVel = efCollision->meshVelocity( dt ) + offset;
            if( addExternalContact( strand1->getGlobalIndex(), edgeIdx1, efCollision->abscissa(), collision ) ){
                ++nInt;
            }
        }

    }
    std::cout << "CTCD " << nCTCD << " ext: " << nInt << std::endl;
}

bool Simulation::addExternalContact( const unsigned strIdx, const unsigned edgeIdx, const Scalar abscissa, const CollidingPair& externalContact )
{
    if( !acceptsCollision( *m_strands[strIdx], edgeIdx, abscissa ) ){
        return false;
    }

    // Discard collisions if their normal is almost parallel to the strand edge and would cause stretching
    const int prevEdge = abscissa == 0. ? edgeIdx - 1 : edgeIdx;
    const Vec3 edge = m_strands[strIdx]->getCurrentTangent( prevEdge );
    if( edge.dot( externalContact.m_normal ) > ALMOST_PARALLEL_COS ){
        return false;
    }

#pragma omp critical
    {
        m_externalContacts[strIdx].push_back( externalContact );
    }

    return true;
}

void Simulation::computeCollidingGroups( const CollidingPairs &origMutualCollisions )
{
    if( origMutualCollisions.empty() ){
        return;
    }

    CollidingPairs mutualCollisions;
    mutualCollisions.reserve( origMutualCollisions.size() );

    for( unsigned i = 0; i < origMutualCollisions.size(); ++i )
    {
        const CollidingPair &proxyCol = origMutualCollisions[i];
        const unsigned s1 = proxyCol.objects.first.globalIndex;
        const unsigned s2 = proxyCol.objects.second.globalIndex;

        if( m_steppers[s1]->refusesMutualContacts() || m_steppers[s2]->refusesMutualContacts() )
        {
            CollidingPair copy( proxyCol );
            makeExternalContact( copy, m_steppers[s2]->refusesMutualContacts() );
        }
        else
        {
            mutualCollisions.push_back( proxyCol );
        }
    }
    std::sort( mutualCollisions.begin(), mutualCollisions.end() );

    // For each strand, list all other contacting ones
    std::vector< std::deque<unsigned> > objsGroups( m_strands.size() );
    for( unsigned i = 0; i < mutualCollisions.size(); ++i )
    {
        const CollidingPair &proxyCol = mutualCollisions[i];
        const unsigned s1 = proxyCol.objects.first.globalIndex;
        const unsigned s2 = proxyCol.objects.second.globalIndex;

        objsGroups[s1].push_back( s2 );
        objsGroups[s2].push_back( s1 );
    }

    // Extract connected subgraphs from the global constraint graph
    for( unsigned s1 = 0; s1 < objsGroups.size(); ++s1 )
    {
        if( m_collidingGroupsIdx[s1] != -1 || objsGroups[s1].empty() ){
            continue;
        }

        const unsigned groupIdx = m_collidingGroups.size();

        m_collidingGroups.push_back( CollidingGroup() );
        CollidingGroup &cg = m_collidingGroups.back();

        m_collidingGroupsIdx[s1] = groupIdx;
        cg.first[s1] = 0;

        std::deque<unsigned> toVisit = objsGroups[s1];

        while( toVisit.size() )
        {
            const unsigned s2 = toVisit.front();
            toVisit.pop_front();

            if( m_collidingGroupsIdx[s2] != -1 ){
                continue;
            }

            m_collidingGroupsIdx[s2] = groupIdx;
            cg.first[s2] = 0;

            toVisit.insert( toVisit.end(), objsGroups[s2].begin(), objsGroups[s2].end() ); // add anything this object is touching
        }

    }

    for( unsigned i = 0; i < mutualCollisions.size(); ++i )
    { // assign the collisions to the groups they belong in
        const CollidingPair &mutualCollision = mutualCollisions[i];
        const unsigned s1 = mutualCollision.objects.first.globalIndex;

        m_collidingGroups[ m_collidingGroupsIdx[s1] ].second.push_back( mutualCollision );
    }

    // Index in group needs to be the same that when using std::map iterator
#pragma omp parallel for
    for ( unsigned i = 0; i < m_collidingGroups.size(); ++i )
    {
        unsigned k = 0;
        IndicesMap &indices = m_collidingGroups[i].first;
        for ( IndicesMap::iterator it = indices.begin(); it != indices.end(); ++it )
        {
            it->second = k++; // assign increasing group indices from 0 -> k
        }
    }
    std::cout << "numCollidingGroups: " << m_collidingGroups.size() << std::endl;
}

void Simulation::setupDeformationBasis( CollidingPair &collision ) const
{
    collision.generateTransformationMatrix();
    computeDeformationGradient( collision.objects.first );
    if( collision.objects.second.globalIndex != -1 ){
        computeDeformationGradient( collision.objects.second );
    }
}

void Simulation::computeDeformationGradient( CollidingPair::Object &object ) const
{
    object.defGrad = new DeformationGradient( 3, m_strands[object.globalIndex]->getCurrentDegreesOfFreedom().rows() );
    DeformationGradient& H = *(object.defGrad);
    H.reserve( object.abscissa > 0. ? 6 : 3 );
    for ( unsigned k = 0; k < 3; ++k )
    {
        H.startVec( k );
        const unsigned col = 4 * object.vertex + k;
        H.insertBackByOuterInner( k, col ) = ( 1. - object.abscissa );

        if ( object.abscissa > 0. )
        {
            H.insertBack( k, col + 4 ) = object.abscissa;
        }
    }

    H.finalize();
}

bool Simulation::needsExternalSolve( unsigned strandIdx ) const
{
    return !m_externalContacts[strandIdx].empty();
}

void Simulation::solveCollidingGroup( CollidingGroup& collisionGroup, bool asFailSafe, bool nonLinear )
{
    if ( collisionGroup.first.empty() ) return;

    std::vector<unsigned> globalIds;

    bool accept = false;
    {
        VecXx vels;
        VecXx impulses;
        bogus::MecheFrictionProblem mecheProblem;
        VecXu startDofs;
        VecXu nDofs;
        if( assembleBogusFrictionProblem( collisionGroup, mecheProblem, globalIds, vels, impulses, startDofs, nDofs ) ){
            accept = solveBogusFrictionProblem( mecheProblem, globalIds, asFailSafe, nonLinear, vels, impulses );
            postProcessBogusFrictionProblem( accept, collisionGroup, mecheProblem, globalIds, vels, impulses, startDofs, nDofs );
        }
    }

    bool mustRetry = !accept;
    if( accept ){
        for ( unsigned i = 0; i < globalIds.size(); ++i )
        {
            const unsigned sIdx = globalIds[i];
            if( m_steppers[sIdx]->lastStepWasRejected() )
            {
                mustRetry = true;
                break;
            }
        }        
    }

    if( mustRetry )
    {
        if( globalIds.size() > 1 )
        { // Failed, drop rod-rod collisions, keep external
#pragma omp parallel for
            for ( unsigned i = 0; i < globalIds.size(); ++i )
            {
                const unsigned sIdx = globalIds[i];
                if( !accept || ( m_steppers[sIdx]->lastStepWasRejected() ) ){
                    solveOnlyStrandExternal( globalIds[i], true, m_params.m_useNonLinearAsFailsafe || m_params.m_alwaysUseNonLinear );
                }
            }
        }
        else if( asFailSafe && nonLinear )
        { // H, otherwise just solve unconstrained and ignore contacts
            ImplicitStepper& stepper = *m_steppers[globalIds[0]];
            stepper.resetStep();
            stepper.solveUnconstrained( true );
            stepper.update();
        }
    }
    
/*
    if ( !asFailSafe )
    { // cleanup some memory
#pragma omp parallel for
        for ( unsigned i = 0; i < colPointers.size(); ++i )
        {
            CollidingPair& col = *colPointers[i];
            delete col.objects.first.defGrad;
            col.objects.first.defGrad = NULL;
            if ( col.objects.second.globalIndex != -1 )
            {
                delete col.objects.second.defGrad;
                col.objects.second.defGrad = NULL;
            }
        }
    }
    */
}

void Simulation::solveOnlyStrandExternal( unsigned objectIdx, bool asFailSafe, bool nonLinear )
{ // send an empty group to solve, since External contacts are added separately in assembly
    CollidingGroup collisionGroup;
    collisionGroup.first[objectIdx] = 0;
    solveCollidingGroup( collisionGroup, asFailSafe, nonLinear );
}

