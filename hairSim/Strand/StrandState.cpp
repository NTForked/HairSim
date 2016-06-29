#include "StrandState.hh"
#include "ElasticStrandUtils.hh"
#include "BandMatrix.hh"

#include "../Utils/Distances.hh"
#include "../Utils/TextLog.hh"

#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <Eigen/Sparse>

StrandState::StrandState( const VecXx& dofs, BendingMatrixBase& bendingMatrixBase ) :
        m_numVertices( dofs.getNumVertices() ),
        m_totalEnergy( 0 ),
        m_dofs( dofs ),
        m_edges( m_dofs ),
        m_lengths( m_edges ),
        m_tangents( m_edges, m_lengths ),
        m_referenceFrames1( m_tangents ),
        m_referenceFrames2( m_tangents, m_referenceFrames1 ),
        m_referenceTwists( m_tangents, m_referenceFrames1 ),
        m_twists( m_referenceTwists, m_dofs ),
        m_curvatureBinormals( m_tangents ),
        m_trigThetas( m_dofs ),
        m_materialFrames1( m_trigThetas, m_referenceFrames1, m_referenceFrames2 ),
        m_materialFrames2( m_trigThetas, m_referenceFrames1, m_referenceFrames2 ),
        m_kappas( m_curvatureBinormals, m_materialFrames1, m_materialFrames2 ),
        m_gradKappas( m_lengths, m_tangents, m_curvatureBinormals, m_materialFrame m_materialFrames2, m_kappas ),
        m_gradTwists( m_lengths, m_curvatureBinormals ),
        m_gradTwistsSquared( m_gradTwists ),
        m_hessKappas( m_lengths, m_tangents, m_curvatureBinormals, m_materialFrame m_materialFrames2, m_kappas ),
        m_hessTwists( m_tangents, m_lengths, m_curvatureBinormals ),
        m_thetaHessKappas( m_curvatureBinormals, m_materialFrames1, m_materialFrames2 ),
        m_bendingProducts( bendingMatrixBase, m_gradKappas ),
{}

StrandState::~StrandState()
{}

void StrandState::resizeSelf()
{
    const size_t ndofs = getDegreesOfFreedom().size();
    assert( ndofs % 4 == 3 );
    // dofs are 3 per vertex, one per edge
    assert( ndofs > 3 );
    // minimum two vertices per strand

    m_totalForce.resize( ndofs );
}

void StrandState::freeCachedQuantities()
{
    m_curvatureBinormals.free();
    m_trigThetas.free();
    m_gradKappas.free();
    m_gradTwists.free();
    m_gradTwistsSquared.free();
    m_hessKappas.free();
    m_hessTwists.free();
    m_thetaHessKappas.free();
    m_bendingProducts.free();
}

bool StrandState::hasSmallForces( const Scalar lTwoTol, const Scalar lInfTol ) const
{
    return ( ( m_totalForce.norm() / m_numVertices <= lTwoTol )
          || ( m_totalForce.lpNorm<Eigen::Infinity>() <= lInfTol ) );
}

Vec3 StrandState::closestPoint( const Vec3& x ) const
{
    Scalar mindist = std::numeric_limits<Scalar>::max();
    Vec3 winner;

    for ( int vtx = 0; vtx < m_numVertices - 1; ++vtx )
    {
        Vec3 y = ClosestPtPointSegment( x, getVertex( vtx ), getVertex( vtx + 1 ) );
        Scalar dist = ( y - x ).squaredNorm();
        if ( dist < mindist )
        {
            mindist = dist;
            winner = y;
        }
    }

    return winner;
}

