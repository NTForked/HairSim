#include "ElasticStrand.hh"
#include "LinearSolver.hh"
#include "ElasticStrandUtils.hh"

#include "../Forces/ForceAccumulator.hh"
#include "../Forces/StretchingForce.hh"
#include "../Forces/TwistingForce.hh"
#include "../Forces/BendingForce.hh"
#include "../Forces/GravitationForce.hh"
#include "../Forces/StrandStrandForce.hh"

#include "../Dynamic/StrandDynamics.hh"

#include "../Utils/TextLog.hh"
#include "../Utils/EigenSerialization.hh"

#include "../Dependencies/ReferenceFrames.hh"

#include "../Collision/CollisionUtils.hh"
#include "../Collision/OrientedBoundingBox.hh"

#include <boost/serialization/access.hpp>
#include <boost/serialization/version.hpp>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/set.hpp>

ElasticStrand::ElasticStrand( 
    const VecXx& dofs, 
    const ElasticStrandParameters& parameters,
    DOFScriptingController* controller = NULL, 
    int globalIndex, 
    const Vec3& initRefFrame1 ) :
        m_globalIndex( globalIndex ),
        m_parameters( parameters ),
        m_numVertices( m_parameters.getNumVertices() ),
        m_currentState( new StrandState( dofs, m_parameters.getBendingMatrixBase() ) ),
        m_futureState( new StrandState( dofs, m_parameters.getBendingMatrixBase() ) ),
        m_dynamics( NULL ), 
        m_requiresExactJacobian( false ),
        m_activelySimulated( true ),
{
    m_totalJacobian = new JacobianMatrixType;
    m_totalJacobian->resize( static_cast<IndexType>( dofs.size() ), static_cast<IndexType>( dofs.size() ) );

    if( controller == NULL ) controller = new DOFScriptingController();
    createDynamics();
    m_dynamics->setScriptingController( controller );
    
    resizeInternals();
    freezeRestShape( 0, m_numEdges ); // for now the rest shape is the shape in which the strand is created, unless modified later on.
}

ElasticStrand::~ElasticStrand()
{
    // Statics/Dynamic should be deleted before StrandStates
    // ( as they may have stuff stored on them, such as StaticSollisionsInfo )
    delete m_dynamics;
    delete m_currentState;
    delete m_futureState;
}

void ElasticStrand::createDynamics()
{
    if ( !m_dynamics )
    {
        m_dynamics = new StrandDynamics( *this );
    }
}

// To be called on creation
void ElasticStrand::resizeInternals()
{
    m_currentState->resizeSelf();
    m_futureState->resizeSelf();

    const IndexType ndofs = getCurrentDegreesOfFreedom().size();
    m_numVertices = ( ndofs + 1 ) / 4;
    m_numEdges = m_numVertices - 1; // This assumes open topology; allowing for loops would start here...

    m_dynamics->resizeSelf();

    m_restLengths.resize( m_numEdges );
    m_restKappas.resize( m_numEdges );
    m_restTwists.resize( m_numEdges );
    m_vertexMasses.resize( m_numVertices );
    m_VoronoiLengths.resize( m_numVertices );
    m_invVoronoiLengths.resize( m_numVertices );
}

void ElasticStrand::invalidatePhysics()
{
    if ( m_dynamics )
        m_dynamics->invalidatePhysics();
}

void ElasticStrand::invalidateFuturePhysics()
{
    if ( m_dynamics )
        m_dynamics->invalidateFuturePhysics();
}

// Take the current geometry as rest shape
void ElasticStrand::freezeRestShape( unsigned begin, unsigned end, Scalar damping )
{
    // Fix rest lengths
    for ( IndexType vtx = begin; vtx < end; ++vtx )
        m_restLengths[vtx] = ( 1. - damping ) * m_currentState->m_lengths[vtx]
                + damping * m_restLengths[vtx];
    updateEverythingThatDependsOnRestLengths();

    for ( IndexType vtx = begin; vtx < end; ++vtx )
    {
        m_restKappas[vtx] = ( 1. - damping ) * m_currentState->m_kappas[vtx]
                + damping * m_restKappas[vtx];
        m_restTwists[vtx] = ( 1. - damping ) * m_currentState->m_twists[vtx]
                + damping * m_restTwists[vtx];
    }

    invalidatePhysics();
}

// Set rest shape from dofs
void ElasticStrand::setRestShape( const VecXx &dofs, unsigned begin, unsigned end, Scalar damping )
{
    VecXx backup = getCurrentDegreesOfFreedom();
    m_currentState->setDegreesOfFreedom( dofs );

    const Vec3 &initRefFrame = m_currentState->getReferenceFrame1( 0 );
    m_currentState->m_referenceFrames1.storeInitialFrames( initRefFrame );

    freezeRestShape( begin, end, damping );

    for ( IndexType vtx = begin + 1; vtx < end; ++vtx )
    {
        m_restTwists[vtx] = m_restTwists[vtx - 1]
                + clamp2Pi( m_restTwists[vtx] - m_restTwists[vtx - 1] );
    }

    m_currentState->setDegreesOfFreedom( backup );

    invalidateCurrentGeometry();
    invalidatePhysics();
}

void ElasticStrand::updateEverythingThatDependsOnRestLengths()
{
    // Total rest length
    m_totalRestLength = 0.0;
    for ( IndexType vtx = 0; vtx < m_numEdges; ++vtx )
        m_totalRestLength += m_restLengths[vtx];

    // Compute Voronoi lengths
    m_VoronoiLengths[0] = 0.5 * m_restLengths[0];
    for ( IndexType vtx = 1; vtx < m_numEdges; ++vtx )
        m_VoronoiLengths[vtx] = 0.5 * ( m_restLengths[vtx - 1] + m_restLengths[vtx] );
    m_VoronoiLengths[m_numEdges] = 0.5 * m_restLengths[m_numVertices - 2];

    // Compute masses and inverse of Voronoi lengths
    for ( IndexType vtx = 0; vtx < m_numVertices; ++vtx )
    {
        m_vertexMasses[vtx] = m_parameters.getDensity() * m_VoronoiLengths[vtx] * M_PI
                * m_parameters.getRadius( vtx ) * m_parameters.getRadius( vtx );

        m_invVoronoiLengths[vtx] = 1.0 / m_VoronoiLengths[vtx];
    }

    invalidatePhysics();
}

void ElasticStrand::setEdgeRestLength( const IndexType vtx, const Scalar newrestlength )
{
    assert( vtx < m_numEdges );

    m_restLengths[vtx] = newrestlength;

    if ( 0 == vtx )
    {
        // If we change the rest length of a fixed edge, move the vertices
        setVertex( 1, getVertex( 0 ) + newrestlength * getEdgeVector( 0 ).normalized() );
    }

}

void ElasticStrand::setEdgesRestLength( const Scalar newRestLength )
{
    for ( IndexType vtx = 0; vtx < m_numEdges; ++vtx )
        setEdgeRestLength( vtx, newRestLength );
    updateEverythingThatDependsOnRestLengths();

    invalidatePhysics();
}

void ElasticStrand::setRadius( const Scalar radius_a )
{
    m_parameters.setRadius( radius_a );

    for ( IndexType vtx = 0; vtx < m_numVertices; ++vtx )
    {
        m_vertexMasses[vtx] = m_parameters.getDensity() * m_VoronoiLengths[vtx] * M_PI
                * m_parameters.getRadius( vtx ) * m_parameters.getRadius( vtx );
    }

    invalidatePhysics();
}

void ElasticStrand::setStiffness( const Scalar youngs )
{
    m_parameters.setYoungsModulus( youngs );

    invalidatePhysics();
}

void ElasticStrand::ackParametersChanged()
{
    updateEverythingThatDependsOnRestLengths();
}

void ElasticStrand::setParameters( const ElasticStrandParameters &parameters )
{
    m_parameters = parameters;

    ackParametersChanged();
}

void ElasticStrand::setParameters( double i_radiusA, double i_rootRM,
        double i_tipRM, double i_youngsModulus, double i_shearModulus, double i_density,
        double i_viscosity, double i_airDrag )
{
    m_parameters.setRadius( i_radiusA );
    m_parameters.setYoungsModulus( i_youngsModulus );
    m_parameters.setShearModulus( i_shearModulus );
    m_parameters.setViscosity( i_viscosity );
    m_parameters.setDensity( i_density );
    m_parameters.setAirDrag( i_airDrag );
    m_parameters.setRootRadiusMultiplier( i_rootRM );
    m_parameters.setTipRadiusMultiplier( i_tipRM );

    ackParametersChanged();
}

void ElasticStrand::setFutureDegreesOfFreedom( const VecXx& dof )
{
    getFutureState().setDegreesOfFreedom( dof );

    invalidateFuturePhysics();
}

void ElasticStrand::setCurrentDegreesOfFreedom( const VecXx& dof )
{

    getCurrentState().setDegreesOfFreedom( dof );

    invalidateCurrentGeometry();
}

Scalar ElasticStrand::getUnsignedAngleToMajorRadius( int vtx, const Vec3& vec ) const
{
   if ( vtx + 1 == m_numVertices )
        --vtx;

    const Vec3& edge = getEdgeVector( vtx ).normalized();
    const Vec3& orth = ( vec - vec.dot( edge ) * edge );
    const Scalar north = orth.norm();
    if ( isSmall( north ) )
        return 0.;
    return std::acos( clamp( getMaterialFrame2( vtx ).dot( orth / north ), -1., 1. ) );
}

Scalar ElasticStrand::getSignedAngleToMajorRadius( int vtx, const Vec3& vec ) const
{
    if ( vtx + 1 == m_numVertices )
        --vtx;

    const Vec3& edge = getEdgeVector( vtx ).normalized();
    const Vec3& orth = ( vec - vec.dot( edge ) * edge );

    const Scalar cosa = getMaterialFrame2( vtx ).dot( orth );
    const Scalar sina = -getMaterialFrame1( vtx ).dot( orth );
    return std::atan2( sina, cosa );
}

// Add ForceT's force on theta only to the VecXx
template<typename ForceT>
void ElasticStrand::accumulateEFThetaOnly( Scalar& thetaE, VecXx& thetaF,
        StrandState& geometry ) const
{
    typename ForceT::LocalThetaForceType localF;

    for ( IndexType vtx = ForceT::s_first; vtx < m_numVertices - ForceT::s_last; ++vtx )
    {
        thetaE += ForceT::localEnergy( *this, geometry, vtx );

        ForceT::computeLocal( localF, *this, geometry, vtx );
        ForceT::addInPosition( thetaF, vtx, localF );
    }
}

// Add ForceT's Jacobian on theta only to the tridiagonal matrix
template<typename ForceT>
void ElasticStrand::accumulateJThetaOnly( TriDiagonalMatrixType& thetaJ,
        StrandState& geometry ) const
{
    typename ForceT::LocalThetaJacobianType localJ;

    for ( IndexType vtx = ForceT::s_first; vtx < m_numVertices - ForceT::s_last; ++vtx )
    {
        ForceT::computeLocal( localJ, *this, geometry, vtx );
        ForceT::addInPosition( thetaJ, vtx, localJ );
    }
}

Scalar ElasticStrand::getCurrentTotalLength() const
{
    Scalar totalLength = 0.0;
    for ( IndexType vtx = 0; vtx < m_numEdges; ++vtx )
    {
        totalLength += m_currentState->m_lengths[vtx];
    }

    return totalLength;
}

Scalar ElasticStrand::getFutureTotalLength() const
{
    Scalar totalLength = 0.0;
    for ( IndexType vtx = 0; vtx < m_numVertices - 1; ++vtx )
    {
        totalLength += m_futureState->m_lengths[vtx];
    }

    return totalLength;
}

void ElasticStrand::findBendElasticLimitMax( Scalar& bendElasticLimit ) const
{
    //find the max BendElastic limit per strand
    for ( IndexType vtx = 1; vtx < m_numEdges; ++vtx )
    {
        for ( int i = 0; i < 2; ++i )
        {
            bendElasticLimit = std::max( bendElasticLimit,
                    fabs( Scalar( m_currentState->m_kappas[vtx][i] - m_restKappas[vtx][i] ) ) );
        }
    }
}

void ElasticStrand::applyPlasticDeformation( Scalar stretchElasticLimit, Scalar bendElasticLimit,
        Scalar twistElasticLimit )
{
    // Rest lengths plastic deformation and everything that depend on them: total rest length, Voronoi lengths, masses...
    for ( int vtx = 0; vtx < m_numEdges; ++vtx )
        if ( m_currentState->m_lengths[vtx] > m_restLengths[vtx] + stretchElasticLimit )
            m_restLengths[vtx] = m_currentState->m_lengths[vtx] - stretchElasticLimit;
    updateEverythingThatDependsOnRestLengths();

    for ( IndexType vtx = 0; vtx < m_numEdges; ++vtx )
    {
        for ( int i = 0; i < 2; ++i )
        {
            if ( m_currentState->m_kappas[vtx][i] > m_restKappas[vtx][i] + bendElasticLimit )
            {
                m_restKappas[vtx][i] = m_currentState->m_kappas[vtx][i] - bendElasticLimit;
            }
            if ( m_currentState->m_kappas[vtx][i] < m_restKappas[vtx][i] - bendElasticLimit )
            {
                m_restKappas[vtx][i] = m_currentState->m_kappas[vtx][i] + bendElasticLimit;
            }
        }

        if ( m_currentState->m_twists[vtx] > m_restTwists[vtx] + twistElasticLimit )
        {
            m_restTwists[vtx] = m_currentState->m_twists[vtx] - twistElasticLimit;
        }
        if ( m_currentState->m_twists[vtx] < m_restTwists[vtx] - twistElasticLimit )
        {
            m_restTwists[vtx] = m_currentState->m_twists[vtx] + twistElasticLimit;
        }
    }

    invalidatePhysics();
}

void ElasticStrand::filterFutureGeometryByRestLength( const double epsilon, bool allowCompression )
{

    Vec3 xaN = m_futureState->getVertex( 0 );
    Vec3 xaP = xaN;

    for ( int i = 1; i < m_numVertices; ++i )
    {
        const Vec3& xbN = m_futureState->getVertex( i );

        const Scalar lP = getEdgeRestLength( i - 1 );

        const Scalar lN = ( xbN - xaP ).norm();
        Scalar ratio = 1.;

        if ( !isSmall( lN ) )
        {
            if ( lN > lP + epsilon )
            {
                ratio = ( lP + epsilon ) / lN;
            }
            else if ( !allowCompression && ( lN < lP - epsilon ) )
            {
                ratio = ( lP - epsilon ) / lN;
            }
        }

        // compute and store revised delta
        const Vec3 xbNrev = xaN + ( xbN - xaP ) * ratio;

        m_futureState->setVertex( i, xbNrev );

        xaP = xbN;
        xaN = xbNrev;
    }

    invalidateFuturePhysics();
}

Scalar ElasticStrand::getCurvilinearAbscissa( int vtx, Scalar localAbscissa ) const
{
    Scalar s = 0;

    for( auto i = 0; i < vtx; ++i )
    {
        s += getEdgeRestLength( i );
    }
    if( vtx + 1 < m_numVertices )
    {
        s += localAbscissa * getEdgeRestLength( vtx );
    }

    return s;
}

void ElasticStrand::getLocalAbscissa( const Scalar curvilinearAbscissa, int &vtx,
        Scalar &localAbscissa ) const
{
    for ( localAbscissa = curvilinearAbscissa, vtx = 0;
            vtx + 1 < m_numVertices && m_restLengths[vtx] < localAbscissa; localAbscissa -=
                    m_restLengths[vtx++] )
        ;

    if ( vtx + 1 == m_numVertices )
    {
        localAbscissa = 0;
    }
    else
    {
        localAbscissa /= m_restLengths[vtx];
    }
}

void ElasticStrand::getSegment( unsigned elementID, Vec3 &start, Vec3 &end ) const
{
    assert( elementID + 1 < getNumVertices() );
    start = getVertex( elementID );
    end = getVertex( elementID + 1 );

}

std::ostream& operator<<( std::ostream& os, const ElasticStrand& strand )
{
    const VecXx& dofs = strand.getCurrentDegreesOfFreedom();
    os << '{';
    for ( int i = 0; i < strand.m_numEdges; i++ )
    {
        os << '{' << dofs[4 * i] << ", " << dofs[4 * i + 1] << ", " << dofs[4 * i + 2] << "}, ";
        // os << strand.m_currentState->m_degreesOfFreedom[4 * i + 3] << ', ';
    }
    os << '{' << dofs[4 * ( strand.m_numEdges )] << ", " << dofs[4 * ( strand.m_numEdges ) + 1]
            << ", " << dofs[4 * ( strand.m_numEdges ) + 2] << '}';
    os << '}';

    return os;
}

/**
 * This will be called as a safeguard if relaxThetas fails: put all the non-fixed thetas to zero. If someone has a better idea...
 */
void ElasticStrand::resetThetas()
{
    m_futureState->setThetas( VecXx( m_numEdges ) );
}

bool ElasticStrand::serializeTo( std::ostream & os ) const
{
    std::cerr << "ElasticStrand::serializeTo not supported temporarily" << std::endl;
    return false;
}

bool ElasticStrand::deserializeFrom( std::istream & is )
{
    std::cerr << "ElasticStrand::deserializeFrom not supported temporarily" << std::endl;
    return false;
}

using namespace std;
void printVec3Array( Vec3Array vec )
{
    cout.precision( std::numeric_limits<double>::digits10 + 2);
    for (unsigned i = 0; i < vec.size(); ++i ){
        cout << vec[i].x() << " " << vec[i].y() << " " << vec[i].z() << " ";
    }
    cout << endl;   
}
