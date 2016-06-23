#ifndef STRANDSTATE_H_
#define STRANDSTATE_H_

#include "Definitions.hh"
#include "BandMatrixFwd.hh"
#include "../Dependencies/Twists.hh"
#include "../Dependencies/BendingProducts.hh"

#include <tr1/memory>

/**
 * This class contains all the strand's variable data (i.e. changed by simulation).
 */
class StrandState
{
public:
    explicit StrandState( const VecXx& dofs, BendingMatrixBase& bendingMatrixBase );

    ~StrandState();

    const VecXx& getDegreesOfFreedom() const
    {
        return m_dofs.get();
    }

    void setDegreesOfFreedom( const VecXx& dof )
    {
        m_dofs.set( dof );
    }

    const Eigen::Map<const VecXx, Eigen::Unaligned, Eigen::InnerStride<4> > getThetas() const
    {
        return m_dofs.getThetas();
    }

    void setThetas( const VecXx& thetas, int numberOfFixedThetas = 0 )
    {
        m_dofs.setThetas( thetas, numberOfFixedThetas );
    }

    const Vec3x getVertex( const IndexType vtx ) const
    {
        assert( vtx < m_numVertices );
        return m_dofs.getVertex( vtx );
    }

    //! Gets a vertex somewhere on an edge
    const Vec3x getVertex( const IndexType vtx, Scalar localAbscissa ) const
    {
        const Vec3x v0 = getVertex( vtx );
        if ( localAbscissa > 0. )
        {
            return ( 1 - localAbscissa ) * v0  +  localAbscissa * getVertex( vtx + 1 );
        }
        return v0;
    }

    Vec3x getVertex( const IndexType vtx, const Mat4x &transform ) const
    {
        assert( vtx < m_numVertices );
        Vec4x v;
        v.head<3>() = m_dofs.getVertex( vtx );
        v[3] = 1.;
        return Vec4x( transform * v ).head<3>();
    }

    //! Gets a vertex somewhere on an edge
    const Vec3x getVertex( const IndexType vtx, Scalar localAbscissa, const Mat4x &transform ) const
    {
        assert( vtx < m_numVertices );

        Vec4x v;
        v.head<3>() = getVertex( vtx, localAbscissa );
        v[3] = 1.;
        return Vec4x( transform * v ).head<3>();
    }

    Scalar getTheta( const IndexType vtx ) const
    {
        assert( vtx < m_numVertices - 1 );
        return m_dofs.getTheta( vtx );
    }

    void setTheta( const IndexType vtx, const Scalar newTheta )
    {
        assert( vtx < m_numVertices - 1 );
        m_dofs.setTheta( vtx, newTheta );
    }

    const Vec3x getEdgeVector( const IndexType vtx ) const
    {
        assert( vtx < m_numVertices - 1 );
        return getVertex( vtx + 1 ) - getVertex( vtx );
    }

    inline int numVertices() const
    {
        return m_numVertices;
    }

    void storeInitialFrames( const Vec3x& initRefFrame1 = Vec3x() )
    {
        m_referenceFrames1.storeInitialFrames( initRefFrame1 );
    }

    void setVertex( const IndexType vtx, const Vec3x& coordinates )
    {
        assert( vtx < m_numVertices );
        m_dofs.setVertex( vtx, coordinates );
    }

    const Vec3x getReferenceFrame1( const IndexType vtx ) const
    {
        assert( vtx < m_numVertices - 1 );
        return m_referenceFrames1[vtx];
    }

    void setReferenceFrame1( const IndexType vtx, const Vec3x& vec )
    {
        assert( vtx < m_numVertices - 1 );
        m_referenceFrames1.set( vtx, vec );
    }

    const Vec3x getReferenceFrame2( const IndexType vtx ) const
    {
        assert( vtx < m_numVertices - 1 );
        return m_referenceFrames2[vtx];
    }

    void setReferenceFrame2( const IndexType vtx, const Vec3x& vec )
    {
        assert( vtx < m_numVertices - 1 );
        m_referenceFrames2.set( vtx, vec );
    }

    const Vec3x getMaterialFrame1( const IndexType vtx ) const
    {
        assert( vtx < m_numVertices - 1 );
        return m_materialFrames1[vtx];
    }

    void setMaterialFrame1( const IndexType vtx, const Vec3x& vec )
    {
        assert( vtx < m_numVertices - 1 );
        m_materialFrames1.set( vtx, vec );
    }

    const Vec3x getMaterialFrame2( const IndexType vtx ) const
    {
        assert( vtx < m_numVertices - 1 );
        return m_materialFrames2[vtx];
    }

    void setMaterialFrame2( const IndexType vtx, const Vec3x& vec )
    {
        assert( vtx < m_numVertices - 1 );
        m_materialFrames2.set( vtx, vec );
    }

    void resizeSelf();
    void freeCachedQuantities();
    bool hasSmallForces( const Scalar lTwoTol, const Scalar lInfTol ) const;
    void computeDeformationGradient( const unsigned edge, const Scalar alpha,
            const VecXx &velocities, SparseRowMatx* &pH, Vec3x &freeVel ) const;
    Vec3x closestPoint( const Vec3x& x ) const;


private:
    // Convenience copy of the original number of vertices (owned by the strand).
    const IndexType m_numVertices;

    // Energy, force
    Scalar m_totalEnergy;
    VecXx m_totalForce;

    ////////////////////////////////////////

    DOFs m_dofs;
    Edges m_edges;
    Lengths m_lengths;
    Tangents m_tangents;
    mutable ReferenceFrames1 m_referenceFrames1;
    mutable ReferenceFrames2 m_referenceFrames2;
    ReferenceTwists m_referenceTwists;
    Twists m_twists;
    CurvatureBinormals m_curvatureBinormals;
    TrigThetas m_trigThetas;
    mutable MaterialFrames<1> m_materialFrames1;
    mutable MaterialFrames<2> m_materialFrames2;
    Kappas m_kappas;
    GradKappas m_gradKappas;
    GradTwists m_gradTwists;
    GradTwistsSquared m_gradTwistsSquared;
    HessKappas m_hessKappas;
    HessTwists m_hessTwists;
    ThetaHessKappas m_thetaHessKappas;
    BendingProducts m_bendingProducts;

};

#endif /* STRANDSTATE_HH_ */
