#include "TwistEdgeHandler.hh"

#include "ElementProxy.hh"
#include "../Core/ElasticStrand.hh"
#include "../Dynamic/StrandDynamicTraits.hh"
#include "../Utils/Distances.hh"
#include "CollisionUtils.hh"

#define CROSS_RESULT_EPSILON 1e-16
#define FRACTIONAL_INFLUENCE_MULTIPLIER 0.25
#define RAD_TO_DEG 57.2957795131

using namespace std;

TwistEdgeHandler::TwistEdgeHandler( const double& thickness ): 
m_frozenCheck(0), m_frozenScene(false), m_thickness( 2.0 * thickness )
{
    qs_prime.first = -1;
    qe_prime.first = -1;
}

TwistEdgeHandler::~TwistEdgeHandler()
{}


void TwistEdgeHandler::getEdgeVerts( const TwistEdge* edge, const bool& startOfStep, Vec3x& firstVert, Vec3x& secondVert )
{
	if( edge->isTwistedBand )
	{

        if( startOfStep && qs_prime.first == edge->uniqueID ){
            firstVert = qs_prime.second.first;
            secondVert = qs_prime.second.second;
            return;
        }
        else if( !startOfStep && qe_prime.first == edge->uniqueID ){
            firstVert = qe_prime.second.first;
            secondVert = qe_prime.second.second;
            return;
        }


		Vec3x alphaFirstVert, alphaSecondVert;
		getEdgeVerts( edge->parents.first, startOfStep, alphaFirstVert, alphaSecondVert );

		Vec3x betaFirstVert, betaSecondVert;
		getEdgeVerts( edge->parents.second, startOfStep, betaFirstVert, betaSecondVert );

        // cout << "edge: " << edge->uniqueID << " alphaFirstVert: " << alphaFirstVert.norm() << " alphaSecondVert: " << alphaSecondVert.norm() << " betaFirstVert: " << betaFirstVert.norm() << " betaSecondVert: " << betaSecondVert.norm() << endl;

        double alpha, beta;
        ClosestPtSegmentSegment( alphaFirstVert, alphaSecondVert, betaFirstVert, betaSecondVert, alpha, beta, firstVert, secondVert );

        double eps = 1e-6;
        if( alpha < eps ) alpha = eps;
        if( alpha > (1.0 - eps) ) alpha = 1.0 - eps;
        firstVert = (1 - alpha ) * alphaFirstVert + alpha * alphaSecondVert;

        if( beta < eps ) beta = eps;
        if( beta > (1.0 - eps) ) beta = 1.0 - eps;
        secondVert = (1 - beta ) * betaFirstVert + beta * betaSecondVert;         

		return;
	}

    // Otherwise original strandVert, we know its position from the strandInfo:
    assert( edge->m_vertexIndex > -1 );
    firstVert = edge->m_strand.getVertex( edge->m_vertexIndex ); // assumed end of timestep position
    if( startOfStep ) firstVert -= edge->m_strand.dynamics().getDisplacement( edge->m_vertexIndex );
    secondVert = edge->m_strand.getVertex( edge->m_vertexIndex + 1 ); // assumed end of timestep position
    if( startOfStep ) secondVert -= edge->m_strand.dynamics().getDisplacement( edge->m_vertexIndex + 1 );
}

void TwistEdgeHandler::getEdgeVerts( const TwistEdge* edge, const bool& startOfStep, Vec3x& firstVert, Vec3x& secondVert, const VecXx& vels, const double& dt  )
{
        std::cout << "call made" << std::endl;

    if( edge->isTwistedBand )
    {
        Vec3x alphaFirstVert, alphaSecondVert;
        getEdgeVerts( edge->parents.first, startOfStep, alphaFirstVert, alphaSecondVert );

        Vec3x betaFirstVert, betaSecondVert;
        getEdgeVerts( edge->parents.second, startOfStep, betaFirstVert, betaSecondVert );

        double alpha, beta;
        ClosestPtSegmentSegment( alphaFirstVert, alphaSecondVert, betaFirstVert, betaSecondVert, alpha, beta, firstVert, secondVert );

        double eps = 1e-6;
        if( alpha < eps ) alpha = eps;
        if( alpha > (1.0 - eps) ) alpha = 1.0 - eps;
        firstVert = (1 - alpha ) * alphaFirstVert + alpha * alphaSecondVert;

        if( beta < eps ) beta = eps;
        if( beta > (1.0 - eps) ) beta = 1.0 - eps;
        secondVert = (1 - beta ) * betaFirstVert + beta * betaSecondVert;         

        return;
    }

    std::cout << "vels actually gets used" << std::endl;

    // Otherwise original strandVert, we know its position from the strandInfo:
    assert( edge->m_vertexIndex > -1 );
    firstVert = edge->m_strand.getVertex( edge->m_vertexIndex ); // assumed end of timestep position
    if( startOfStep ) firstVert -= ( vels.segment<3>(4 * edge->m_vertexIndex) * dt ); 
    secondVert = edge->m_strand.getVertex( edge->m_vertexIndex + 1 ); // assumed end of timestep position
    if( startOfStep ) secondVert -= ( vels.segment<3>(4 * (edge->m_vertexIndex + 1 ) ) * dt ); 
}


void TwistEdgeHandler::getEdgeAlphas( const TwistEdge* edge, const bool& startOfStep, double& alpha, double& beta )
{
    assert( edge->isTwistedBand );

    Vec3x alphaFirstVert, alphaSecondVert;
    getEdgeVerts( edge->parents.first, startOfStep, alphaFirstVert, alphaSecondVert );

    Vec3x betaFirstVert, betaSecondVert;
    getEdgeVerts( edge->parents.second, startOfStep, betaFirstVert, betaSecondVert );

    Vec3x firstVert, secondVert;
    ClosestPtSegmentSegment( alphaFirstVert, alphaSecondVert, betaFirstVert, betaSecondVert, alpha, beta, firstVert, secondVert );
}

bool TwistEdgeHandler::computeTunnSign( TwistEdge* first, TwistEdge* second, const bool& startOfStep )
{ // signed tet volume: triple product, (ba) cross (ca) dot (da)
    Vec3x a, b, c, d;
    getEdgeVerts( first, startOfStep, a, b );
    getEdgeVerts( second, startOfStep, c, d );
    double crossResult = ((b - a).cross(c - a)).dot( d - a );
    // cout << "["<< first->uniqueID << "-"<< second->uniqueID <<"] crossResult: " << crossResult << " bool: " << (crossResult > CROSS_RESULT_EPSILON) << endl;
    return crossResult > CROSS_RESULT_EPSILON;
}

template <typename T> int sgn(T val)
{ //  returns 1 : pos, -1 : neg, 0 : zero
    return (T(0) < val) - (val < T(0));
}

bool signPos( double number )
{
    return number > 1e-17; // tweak this..
}

void TwistEdgeHandler::updateTwistAngle( TwistEdge* edge, TwistEdge* startPA, TwistEdge* startPB, TwistEdge* endPA, TwistEdge* endPB, const bool& traversal )
{
// cout << "edge: " << edge->uniqueID << " sA: " << startPA->uniqueID << " sB: " << startPB->uniqueID << " eA: " << endPA->uniqueID << " eB: " << endPB->uniqueID << endl;

    Vec3x edgeA, edgeB;
    getEdgeVerts( edge, true, edgeA, edgeB );
    Vec3x planeNormal = edgeB - edgeA;
    planeNormal.normalize();

    Vec3x aS1, aS2; // alpha start first, second
    Vec3x bS1, bS2;
    getEdgeVerts( startPA, true, aS1, aS2 );
    getEdgeVerts( startPB, true, bS1, bS2 );

    Vec3x a = aS2 - aS1; a = a - a.dot(planeNormal) * planeNormal; // project each onto plane
    Vec3x b = bS2 - bS1; b = b - b.dot(planeNormal) * planeNormal;
    double startAngle = atan2( ( a.cross(b) ).dot(planeNormal), a.dot(b) );

    getEdgeVerts( edge, false, edgeA, edgeB );
    planeNormal = edgeB - edgeA;
    planeNormal.normalize();

    Vec3x aF1, aF2; // alpha final first, second
    Vec3x bF1, bF2;
    getEdgeVerts( endPA, false, aF1, aF2 );
    getEdgeVerts( endPB, false, bF1, bF2 );

    a = aF2 - aF1; a = a - a.dot(planeNormal) * planeNormal;
    b = bF2 - bF1; b = b - b.dot(planeNormal) * planeNormal;
    double finalAngle = atan2( ( a.cross(b) ).dot(planeNormal), a.dot(b) );

//     if( edge->currAngle() == 365 )
//     { // When created (first time through...)
//         edge->intersections.top()->originalTwistAngle = finalAngle;
//         edge->intersections.top()->currAngle = finalAngle;

// cout << "first INTERSECTION twist" << endl;
// cout << "intersectionTwists: " << edge->intersectionTwists()  <<  " coplanarTwists: " << edge->coplanarTwists() << " startAngle: " << startAngle * RAD_TO_DEG << " finalAngle: " << finalAngle * RAD_TO_DEG << endl;
//         return;
//     }

    double times[4];
    unsigned num_times;
    getCoplanarityTimes( aS1, aS2, bS1, bS2, aF1, aF2, bF1, bF2, times, NULL, num_times );
    bool intersecting = false;
    if( num_times > 0 ){
        // std::cout << "num_times: " << num_times << std::endl;
        if( times[0] < 1e-12 ) return;

        double alpha, beta;
        Vec3x colAlpha, colBeta;
        Vec3x aC1 = aS1 + times[0] * ( aF1 - aS1 );
        Vec3x aC2 = aS2 + times[0] * ( aF2 - aS2 );
        Vec3x bC1 = bS1 + times[0] * ( bF1 - bS1 );
        Vec3x bC2 = bS2 + times[0] * ( bF2 - bS2 );
        double proxsq = ClosestPtSegmentSegment( aC1, aC2, bC1, bC2, alpha, beta, colAlpha, colBeta );
        // cout << "proxsq " << proxsq << " time: " << times[0] << endl;
        intersecting = proxsq < 1e-12; // -16

        if( !traversal && intersecting )
        {
            if( edge->coplanarTwists() == 0 && edge->intersectionTwists() > 0 ){
                edge->intersections.pop();
                cout << "INTERSECTION twist, popping from STACK " << endl;    
            }
            else{
                edge->intersections.push( new TwistIntersection( finalAngle ) );
                // cout << "INTERSECTION twist, increasing STACK " << endl;
                    edge->intersections.top()->originalTwistAngle = finalAngle;
                    edge->intersections.top()->currAngle = finalAngle;
    // if( startPA->m_strand.getGlobalIndex() / 33 != startPB->m_strand.getGlobalIndex() / 33 ){

                    // cout << "first INTERSECTION twist" << endl;
                    // cout << "intersectionTwists: " << edge->intersectionTwists() <<  " coplanarTwists: " << edge->coplanarTwists() << " startAngle: " << startAngle * RAD_TO_DEG << " finalAngle: " << finalAngle * RAD_TO_DEG << endl;
                    // }
                    return;                
            }
        }
    }
    if( edge->intersectionTwists() == 0 ) return;

    // if( edge->intersectionTwists % 2 == 0 ){
    //     finalAngle = -finalAngle;
    //     if( !intersecting ){
    //         startAngle = -startAngle;
    //     }
    // }

    double deltaAngle = finalAngle - startAngle;

    if( deltaAngle > 3.14159 ){
        cout << "deltaAngle greater than 180: " << deltaAngle * RAD_TO_DEG << endl;
        deltaAngle -= 6.28319;
    }
    else if( deltaAngle < -3.14159 ){
        cout << "deltaAngle less than 180: " << deltaAngle * RAD_TO_DEG << endl;        
        deltaAngle += 6.28319;
    }

// cout << "intersectionTwists: " << edge->intersectionTwists() << " " << edge->intersections.size() << endl;
// cout << "start difference between the two angles : " << ( edge->currAngle() - startAngle ) * RAD_TO_DEG <<endl;

    if( signPos( edge->currAngle() ) != signPos( finalAngle ) && std::abs( edge->currAngle() - finalAngle ) > CROSS_RESULT_EPSILON ){
        cout << "signPos currAngle: " << signPos( edge->currAngle() ) << " finalAngle: " << signPos( finalAngle ) << endl;
        if( signPos(finalAngle) == 1 ){
            edge->intersections.top()->coplanarTwists += 1;
        }
        else if( signPos( edge->currAngle() ) == 1 ){
            edge->intersections.top()->coplanarTwists -= 1;
        }
    }
    else if( ((int) std::abs( edge->currAngle() * RAD_TO_DEG )) / 180 != ((int)std::abs( finalAngle * RAD_TO_DEG )) / 180 ){
        cout << "modulus: " << (int) std::abs(edge->currAngle() * RAD_TO_DEG) << " " << (int) std::abs(finalAngle * RAD_TO_DEG) << endl;
        if( finalAngle > edge->currAngle() ){
            edge->intersections.top()->coplanarTwists += 1;
        }
        else{
            edge->intersections.top()->coplanarTwists -= 1;
        }
    }

    if( intersecting )
    {
        edge->intersections.top()->currAngle = finalAngle;
    }
    else{
        edge->intersections.top()->currAngle += deltaAngle;
    }

// cout << "intersectionTwists: "<< edge->intersectionTwists() << " coplanarTwists: " << edge->coplanarTwists() << " startAngle: " << startAngle * RAD_TO_DEG << " finalAngle: " << finalAngle * RAD_TO_DEG << endl;
// cout << "difference between the two angles: " << (edge->currAngle() - finalAngle) * RAD_TO_DEG <<endl;

//     if( startPA->m_strand.getGlobalIndex() / 33 != startPB->m_strand.getGlobalIndex() / 33 ){
//         std::cout << startPA->m_strand.getGlobalIndex() << " " << startPB->m_strand.getGlobalIndex() << std::endl;
// cout << "intersectionTwists: "<< edge->intersectionTwists() << " coplanarTwists: " << edge->coplanarTwists();
// cout << " originalTwistAngle: " << edge->originalTwistAngle() * RAD_TO_DEG << " currAngle: " << edge->currAngle() * RAD_TO_DEG << endl;        
//     }

}

void TwistEdgeHandler::deleteBand( TwistEdge* edge, std::vector< ElementProxy* >& elementProxies )
{
	// Can only delete leaf nodes
	if( edge->children.size() != 0 ) return;

    // Delete myself from parent's lists:

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

    // Delete myself from tunneled bands list
    std::vector< TwistEdge* >::iterator it = std::find( tunnelingBands.begin(), tunnelingBands.end(), edge );
    if( it != tunnelingBands.end() ){ 
    	tunnelingBands.erase( it );
    }
    else{ 
        cerr << "WARNING, EDGE NOT FOUND TO BE DELETED on tunnelingBands: " << edge->uniqueID << endl;
    }

    // Delete myself from collision detector's proxies
    std::vector< ElementProxy* >::iterator eit = std::find( elementProxies.begin(), elementProxies.end(), edge );
    if( eit != elementProxies.end() ){ 
    	elementProxies.erase( eit );
    }
    else{ 
        cerr << "WARNING, EDGE NOT FOUND TO BE DELETED from elementProxies: " << edge->uniqueID << endl;
    }
}

void TwistEdgeHandler::applyImpulses( ElasticStrand* strand, ImplicitStepper* stepper, bool computeImpulses )
{
    stepper->m_impulseChanges.resize( strand->dynamics().getDisplacements().size() );
    stepper->m_impulseChanges.setZero(); // changes to velocities

    bool impulsesON = true;
    if( !impulsesON || !computeImpulses ){    return;    }

    // std::cout << "applyImpulses " << tunnelingBands.size() << std::endl;

    int sidx = strand->getGlobalIndex();
    double twistScale = 0.0; //5.5;
    double springScale = 10.00;
    double dampingCoeff = 0.20; // 0.35;

    VecXx velocities = stepper->m_strand.dynamics().getDisplacements() / stepper->m_dt;
    strand->rod_data->m_renderer->displaceVec = stepper->m_impulseChanges;

    for( unsigned tb = 0; tb < tunnelingBands.size(); ++tb )
    {
        VecXx prevImpulseChanges = stepper->m_impulseChanges;
        TwistEdge* edge = tunnelingBands[tb];

        if( edge->parents.first->m_strand.getGlobalIndex() != sidx && 
            edge->parents.second->m_strand.getGlobalIndex() != sidx ){
            continue;
        }


        Vec3x edgeA, edgeB;
        // getEdgeVerts( edge, true, edgeA, edgeB ); // should be true if used before step_dynamics, false otherwise
        getEdgeVerts( edge, true, edgeA, edgeB, velocities, stepper->m_dt );// should be true if used before step_dynamics, false otherwise
        Vec3x planeNormal = edgeA - edgeB;

        double x = planeNormal.norm();
        planeNormal.normalize();

        if( edge->intersectionTwists() != 0 ){
            x += 2 * edge->radius;
            x *= 2.5;
            planeNormal *= -1.0;
        }
        if( edge->intersectionTwists() == 0 ) x = (2*edge->radius) - x;
        
        if( x < 0.0 ){
            x = 0.0;
            dampingCoeff = x;
        }

        if( edge->parents.first->m_strand.getGlobalIndex() == sidx )
        {
            int v = edge->parents.first->m_vertexIndex;


            double along1 = dampingCoeff * velocities.segment<3>(4*v).dot( planeNormal );
            double along2 = dampingCoeff * velocities.segment<3>(4 * (v+1) ).dot( planeNormal );
            stepper->m_impulseChanges.segment<3>( 4 * v ) += dampingCoeff * -( velocities.segment<3>(4*v) ) + along1 * planeNormal;
            stepper->m_impulseChanges.segment<3>( 4 * (v + 1) ) += dampingCoeff * -( velocities.segment<3>( 4 * (v+1) ) ) + along2 * planeNormal;

            double spring = springScale * ( x );
            stepper->m_impulseChanges.segment<3>( 4 * v ) += spring * planeNormal;
            stepper->m_impulseChanges.segment<3>( 4 * (v + 1) ) += spring * planeNormal; 
        }

        if( edge->parents.second->m_strand.getGlobalIndex() == sidx )
        {
            int v = edge->parents.second->m_vertexIndex;


            double along1 = dampingCoeff * velocities.segment<3>(4*v).dot( -planeNormal );
            double along2 = dampingCoeff * velocities.segment<3>(4 * (v+1) ).dot( -planeNormal );
            stepper->m_impulseChanges.segment<3>( 4 * v ) += dampingCoeff * -( velocities.segment<3>(4*v) ) + along1 * -planeNormal;
            stepper->m_impulseChanges.segment<3>( 4 * (v + 1) ) += dampingCoeff * -( velocities.segment<3>( 4 * (v+1) ) ) + along2 * -planeNormal;;

            double spring = -springScale * ( x );
            stepper->m_impulseChanges.segment<3>( 4 * v ) += spring * planeNormal;
            stepper->m_impulseChanges.segment<3>( 4 * (v + 1) ) += spring * planeNormal;
        }

        velocities += (stepper->m_impulseChanges - prevImpulseChanges);
        // strand->rod_data->m_renderer->displaceVec += stepper->m_impulseChanges;
        // std::cout << "deltaPenalty: " << (stepper->m_impulseChanges - prevImpulseChanges).transpose() << std::endl;
    }

}

