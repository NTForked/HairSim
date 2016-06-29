#include "TwistEdgeHandler.hh"

#include "ElementProxy.hh"
#include "../Core/ElasticStrand.hh"
#include "../Dynamic/StrandDynamics.hh"
#include "../Utils/Distances.hh"
#include "CollisionUtils.hh"

#define CROSS_RESULT_EPSILON 1e-16
#define FRACTIONAL_INFLUENCE_MULTIPLIER 0.25
#define RAD_TO_DEG 57.2957795131

#define TOO_CLOSE_TO_VERT 1e-6

using namespace std;

TwistEdgeHandler::TwistEdgeHandler(): 
    m_frozenCheck(0), 
    m_frozenScene(false)
{}

TwistEdgeHandler::~TwistEdgeHandler()
{}

void TwistEdgeHandler::getEdgeVerts( const TwistEdge* edge, const bool& startOfStep, Vec3& firstVert, Vec3& secondVert )
{
	if( edge->isTwistedBand )
	{
		Vec3 alphaFirstVert, alphaSecondVert;
		getEdgeVerts( edge->parents.first, startOfStep, alphaFirstVert, alphaSecondVert );

		Vec3 betaFirstVert, betaSecondVert;
		getEdgeVerts( edge->parents.second, startOfStep, betaFirstVert, betaSecondVert );

        double alpha, beta;
        ClosestPtSegmentSegment( alphaFirstVert, alphaSecondVert, betaFirstVert, betaSecondVert, alpha, beta, firstVert, secondVert );

        // clamp points too close to the ends
        double eps = TOO_CLOSE_TO_VERT;
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

    if( startOfStep )
    {
        firstVert = edge->m_strand.getVertex( edge->m_vertexIndex );
        secondVert = edge->m_strand.getVertex( edge->m_vertexIndex + 1 );
    }
    else{
        firstVert = edge->m_strand.getFutureVertex( edge->m_vertexIndex );
        secondVert = edge->m_strand.getFutureVertex( edge->m_vertexIndex + 1 );        
    }
}

void TwistEdgeHandler::getEdgeAlphas( const TwistEdge* edge, const bool& startOfStep, double& alpha, double& beta )
{
    assert( edge->isTwistedBand );

    Vec3 alphaFirstVert, alphaSecondVert;
    getEdgeVerts( edge->parents.first, startOfStep, alphaFirstVert, alphaSecondVert );

    Vec3 betaFirstVert, betaSecondVert;
    getEdgeVerts( edge->parents.second, startOfStep, betaFirstVert, betaSecondVert );

    Vec3 firstVert, secondVert;
    ClosestPtSegmentSegment( alphaFirstVert, alphaSecondVert, betaFirstVert, betaSecondVert, alpha, beta, firstVert, secondVert );
}

bool signPos( double number )
{
    return number > 1e-17; // tweak this..
}

void TwistEdgeHandler::updateTwistAngle( TwistEdge* edge, TwistEdge* startPA, TwistEdge* startPB, TwistEdge* endPA, TwistEdge* endPB, const bool& traversal )
{
    Vec3 edgeA, edgeB;
    getEdgeVerts( edge, true, edgeA, edgeB );
    Vec3 planeNormal = edgeB - edgeA;
    planeNormal.normalize();

    Vec3 aS1, aS2; // alpha start first, second
    Vec3 bS1, bS2;
    getEdgeVerts( startPA, true, aS1, aS2 );
    getEdgeVerts( startPB, true, bS1, bS2 );

    Vec3 a = aS2 - aS1; a = a - a.dot(planeNormal) * planeNormal; // project each onto plane
    Vec3 b = bS2 - bS1; b = b - b.dot(planeNormal) * planeNormal;
    double startAngle = atan2( ( a.cross(b) ).dot(planeNormal), a.dot(b) );

    getEdgeVerts( edge, false, edgeA, edgeB );
    planeNormal = edgeB - edgeA;
    planeNormal.normalize();

    Vec3 aF1, aF2; // alpha final first, second
    Vec3 bF1, bF2;
    getEdgeVerts( endPA, false, aF1, aF2 );
    getEdgeVerts( endPB, false, bF1, bF2 );

    a = aF2 - aF1; a = a - a.dot(planeNormal) * planeNormal;
    b = bF2 - bF1; b = b - b.dot(planeNormal) * planeNormal;
    double finalAngle = atan2( ( a.cross(b) ).dot(planeNormal), a.dot(b) );

    double times[4];
    unsigned num_times;
    getCoplanarityTimes( aS1, aS2, bS1, bS2, aF1, aF2, bF1, bF2, times, NULL, num_times );
    bool intersecting = false;
    if( num_times > 0 ){

        if( times[0] < 1e-12 ) return;

        double alpha, beta;
        Vec3 colAlpha, colBeta;
        Vec3 aC1 = aS1 + times[0] * ( aF1 - aS1 );
        Vec3 aC2 = aS2 + times[0] * ( aF2 - aS2 );
        Vec3 bC1 = bS1 + times[0] * ( bF1 - bS1 );
        Vec3 bC2 = bS2 + times[0] * ( bF2 - bS2 );
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
                edge->intersections.top()->originalTwistAngle = finalAngle;
                edge->intersections.top()->currAngle = finalAngle;
                return;
            }
        }
    }
    if( edge->intersectionTwists() == 0 ) return;

    double deltaAngle = finalAngle - startAngle;

    if( signPos( edge->currAngle() ) != signPos( finalAngle ) && std::abs( edge->currAngle() - finalAngle ) > CROSS_RESULT_EPSILON ){
        if( signPos(finalAngle) == 1 ){
            edge->intersections.top()->coplanarTwists += 1;
        }
        else if( signPos( edge->currAngle() ) == 1 ){
            edge->intersections.top()->coplanarTwists -= 1;
        }
    }
    else if( ((int) std::abs( edge->currAngle() * RAD_TO_DEG )) / 180 != ((int)std::abs( finalAngle * RAD_TO_DEG )) / 180 ){
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

    int sidx = strand->getGlobalIndex();
    double springScale = 5.00;
    double dampingCoeff = 0.20;

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

        Vec3 edgeA, edgeB;
        getEdgeVerts( edge, false, edgeA, edgeB );
        Vec3 planeNormal = edgeA - edgeB;
        double x = planeNormal.norm();
        planeNormal.normalize();

        if( edge->intersectionTwists() != 0 ){
            x += 2 * edge->m_radius;
            x *= 2.5;
            planeNormal *= -1.0;
        }
        if( edge->intersectionTwists() == 0 ){
            x = ( 2 * edge->m_radius ) - x;
        }
        
        if( x < 0.0 ){
            x = 0.0;
            dampingCoeff = x;
        }

        if( edge->parents.first->m_strand.getGlobalIndex() == sidx )
        {
            int v = edge->parents.first->m_vertexIndex;

            double along1 = dampingCoeff * velocities.segment<3>( 4 * v ).dot( planeNormal );
            double along2 = dampingCoeff * velocities.segment<3>( 4 * (v+1) ).dot( planeNormal );
            stepper->m_futureVelocities.segment<3>( 4 * v ) += dampingCoeff * -( velocities.segment<3>(4*v) ) + along1 * planeNormal;
            stepper->m_futureVelocities.segment<3>( 4 * (v + 1) ) += dampingCoeff * -( velocities.segment<3>( 4 * (v+1) ) ) + along2 * planeNormal;

            double spring = springScale * ( x );
            stepper->m_futureVelocities.segment<3>( 4 * v ) += spring * planeNormal;
            stepper->m_futureVelocities.segment<3>( 4 * (v + 1) ) += spring * planeNormal; 
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

