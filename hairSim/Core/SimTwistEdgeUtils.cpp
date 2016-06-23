void StrandImplicitManager::deleteInvertedProxies()
{
    
    if( penaltyAfter )
    {
        for ( std::vector<ElasticStrand*>::size_type i = 0; i < m_strands.size(); i++ )
        {
            m_collisionDetector->m_proxyHistory->applyImpulses( m_strands[i], m_steppers[i], true ); // gather impulse changes
            m_steppers[i]->newVelocities() += m_steppers[i]->m_impulseChanges; // update velocities
            m_steppers[i]->update( true );
            m_steppers[i]->finalize();
        }
    }


    // grab just the tunneled bands, none of the originals
    std::vector< TwistEdge* > tunneledBands = m_collisionDetector->m_proxyHistory->tunnelingBands;
    for( unsigned t = 0; t < tunneledBands.size(); ++t )
    {

        TwistEdge* edge = tunneledBands[t];

        if( penaltyOnce ){
            m_collisionDetector->m_proxyHistory->deleteBand( edge, m_collisionDetector->m_elementProxies );                
            continue;
        }

        // m_collisionDetector->m_proxyHistory->updateTwistAngle( edge );
        if( edge->traversed ){
            edge->traversed = false;
        }
        else{
            m_collisionDetector->m_proxyHistory->updateTwistAngle( edge, edge->parents.first, edge->parents.second, edge->parents.first, edge->parents.second, false );
        }

        if( edge->intersectionTwists() == 0 )
        {
    // cout << "\nattempting to DELETE: " << edge->uniqueID << endl;

            //check thickness boundary...
            Vec3x edgeA, edgeB;
            m_collisionDetector->m_proxyHistory->getEdgeVerts( edge, false, edgeA, edgeB );

            // if( (edgeB - edgeA).norm() >= 2 * edge->radius ){
            if( (edgeB - edgeA).norm() >= (2.0 * edge->radius) ){
                // cout << "\nDELETING TUNN: " << edge->uniqueID << endl;
                m_collisionDetector->m_proxyHistory->deleteBand( edge, m_collisionDetector->m_elementProxies );                
            }
        }
    }  
}

void StrandImplicitManager::traversalCheck()
{
    // TwistEdgeHandler* teh = m_collisionDetector->m_proxyHistory;
    std::vector< TwistEdge* > tunneledBands = m_collisionDetector->m_proxyHistory->tunnelingBands;

    for( unsigned b = 0; b < tunneledBands.size(); ++b ){

        TwistEdge* edge = tunneledBands[b];

        if( edge->parents.first->isTwistedBand && edge->parents.second->isTwistedBand ) continue;

        bool alphaPair = false;
        bool betaPair = false;
        bool alphaPrev = false;
        bool betaPrev = false;

        Scalar alpha, beta;
        m_collisionDetector->m_proxyHistory->getEdgeAlphas( edge, false, alpha, beta );

        if( !edge->parents.first->isTwistedBand ){
            if( alpha < 0.0001 && edge->parents.first->prev != NULL ){
                alphaPair = true;
                alphaPrev = true;
            }
            else if( alpha > 0.9999 && edge->parents.first->next != NULL ){
                alphaPair = true;
            }
        }
        if( !edge->parents.second->isTwistedBand ){
            if( beta < 0.0001 && edge->parents.second->prev != NULL ){
                betaPair = true;
                betaPrev = true;
            }
            else if( beta > 0.9999 && edge->parents.second->next != NULL ){
                betaPair = true;
            }
        }

        // attempting to change traversal rules... to create symmetry on reverse-scripting
        if( alphaPair && !betaPair )
        {
            if( beta < 0.5 && edge->parents.second->prev != NULL ){
                betaPrev = true;
                betaPair = true;
            }
            if( beta > 0.5 && edge->parents.second->next != NULL ){
                betaPrev = false;
                betaPair = true;
            }
        }
        else if( !alphaPair && betaPair )
        {
            if( alpha < 0.5 && edge->parents.first->prev != NULL ){
                alphaPrev = true;
                alphaPair = true;
            }
            if( alpha > 0.5 && edge->parents.first->next != NULL ){
                alphaPrev = false;
                alphaPair = true;
            }
        }

        // not straddling either parent
        if( !alphaPair && !betaPair ) continue;

        // cout << "edge alpha: " << alpha << " beta: " << beta << endl;

        int alphaMax = alphaPair ? 2 : 1;
        int betaMax = betaPair ? 2 : 1;

        TwistEdge* alphaCheck = alphaPrev ? edge->parents.first->prev : edge->parents.first;
        TwistEdge* betaCheck;
        TwistEdge* minAlpha;
        TwistEdge* minBeta;
        double minDistSq = std::numeric_limits<Scalar>::infinity();
        Vec3x bestFirst, bestSecond;

        for( int a = 0; a < alphaMax; ++a ){

            betaCheck = betaPrev ? edge->parents.second->prev : edge->parents.second;
            for( int b = 0; b < betaMax; ++b ){

                Vec3x alphaFirstVert, alphaSecondVert;
                m_collisionDetector->m_proxyHistory->getEdgeVerts( alphaCheck, false, alphaFirstVert, alphaSecondVert );

                Vec3x betaFirstVert, betaSecondVert;
                m_collisionDetector->m_proxyHistory->getEdgeVerts( betaCheck, false, betaFirstVert, betaSecondVert );

                Vec3x firstVert, secondVert;
                double distsq = ClosestPtSegmentSegment( alphaFirstVert, alphaSecondVert, betaFirstVert, betaSecondVert, alpha, beta, firstVert, secondVert );                
                // cout <<" testing: " << alphaCheck->uniqueID << " vs: " << betaCheck->uniqueID << " result: " <<  distsq << endl;

                if( distsq < minDistSq ){

                    double eps = 1e-6;
                    if( alpha < eps ) alpha = eps;
                    if( alpha > (1.0 - eps) ) alpha = 1.0 - eps;
                    firstVert = (1 - alpha ) * alphaFirstVert + alpha * alphaSecondVert;

                    if( beta < eps ) beta = eps;
                    if( beta > (1.0 - eps) ) beta = 1.0 - eps;
                    secondVert = (1 - beta ) * betaFirstVert + beta * betaSecondVert; 

                    minDistSq = distsq;
                    minAlpha = alphaCheck;
                    minBeta = betaCheck;
                    bestFirst = firstVert;
                    bestSecond = secondVert;
                }

                betaCheck = betaCheck->next;
            }
            alphaCheck = alphaCheck->next;
        }

        if( minAlpha != edge->parents.first || minBeta != edge->parents.second )
        { // TRAVERSE

            // cout << "TRAVERSE" << endl;
            // NOT combining for now, just traversing if free of 'similar' constraint on desired edge:
            bool blocked = false;
            for( unsigned c = 0; c < minAlpha->children.size(); ++c ){
                if( minAlpha->children[c].first == minBeta->uniqueID ){
                    blocked = true;
                    // cout << "BLOCKED" << endl;
                    // combine currAngle, coplanarTwists of top intersections, and combine the rest of the intersection stacks into one
                    // or, alternatively, only combine when number of intersections matches? ..../
                    break;
                }
            }
            if( !blocked ){

                // original start
                Vec3x firstVert, secondVert;
                m_collisionDetector->m_proxyHistory->getEdgeVerts( edge, true, firstVert, secondVert );
                m_collisionDetector->m_proxyHistory->qs_prime.first = edge->uniqueID;
                m_collisionDetector->m_proxyHistory->qs_prime.second.first = firstVert;
                m_collisionDetector->m_proxyHistory->qs_prime.second.second = secondVert;

                // new end position
                m_collisionDetector->m_proxyHistory->qe_prime.first = edge->uniqueID;
                m_collisionDetector->m_proxyHistory->qe_prime.second.first = bestFirst;
                m_collisionDetector->m_proxyHistory->qe_prime.second.second = bestSecond;

                edge->frozenID = m_collisionDetector->m_proxyHistory->m_frozenCheck;

                edge->flagged = true;
                edge->parents.first->flagged = true;
                edge->parents.second->flagged = true;
                minAlpha->flagged = true;
                minBeta->flagged = true;

                m_collisionDetector->m_proxyHistory->updateTwistAngle( edge, edge->parents.first, edge->parents.second, minAlpha, minBeta, true );
                edge->traversed = true;

                setupContinuousTimeCollisions();

                edge->flagged = false;
                edge->parents.first->flagged = false;
                edge->parents.second->flagged = false;
                minAlpha->flagged = false;
                minBeta->flagged = false;

                std::pair< int, TwistEdge* > childEdge( edge->parents.second->uniqueID, edge );

                std::vector< std::pair< int, TwistEdge* > >::iterator pit = std::find( edge->parents.first->children.begin(), edge->parents.first->children.end(), childEdge );
                if( pit != edge->parents.first->children.end() ){ 
                    edge->parents.first->children.erase( pit );
                }
                else{ 
                    cerr << "WARNING, EDGE NOT FOUND TO BE DELETED on parent: " << edge->parents.first->uniqueID << endl;
                }

                childEdge.first = edge->parents.first->uniqueID;
                pit = std::find( edge->parents.second->children.begin(), edge->parents.second->children.end(), childEdge );
                if( pit != edge->parents.second->children.end() ){ 
                    edge->parents.second->children.erase( pit );
                }
                else{ 
                    cerr << "WARNING, EDGE NOT FOUND TO BE DELETED on parent: " << edge->parents.second->uniqueID << endl;
                }


                minAlpha->children.push_back( std::pair< int, TwistEdge* >( minBeta->uniqueID, edge ) );
                minBeta->children.push_back( std::pair< int, TwistEdge* >( minAlpha->uniqueID, edge ) );
               
                m_collisionDetector->m_proxyHistory->qs_prime.first = -1;
                m_collisionDetector->m_proxyHistory->qe_prime.first = -1;

                edge->parents.first = minAlpha;
                edge->parents.second = minBeta;

                // cout << "finished" << endl;
            }
        }
    }   
}
