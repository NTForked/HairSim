#ifndef DEGREESOFFREEDOM_HH_
#define DEGREESOFFREEDOM_HH_

#include "DependencyNode.hh"

/**
 * Unit: cm for position dofs, no dimension for theta
 */
class DOFs: public DependencyNode<VecXx>
{
public:
    DOFs( const VecXx& dofValues ) :
            DependencyNode<VecXx>( dofValues )
    {
        assert( dofValues.size() % 4 == 3 );
        m_numEdges = dofValues.size() / 4;
        setClean();
    }

    // As this class is meant to be the root of the dependency tree, we provide a const getter
    // This means that DOFs are never "computed" but updated by the outer loop algorithm
    const VecXx& get() const
    {
        return m_value;
    }

    VecXx& get() 
    {
        return m_value;
    }

    Vec3x getVertex( IndexType vtx ) const
    {
        assert( vtx < (m_numEdges + 1) );

        return get().segment<3>( 4 * vtx );
    }

    void setVertex( IndexType vtx, const Vec3x& point )
    {
        m_value.segment<3>( 4 * vtx ) = point;
        setDependentsDirty();
    }

    // Accessors to the theta degrees of freedom
    const Eigen::Map<const VecXx, Eigen::Unaligned, Eigen::InnerStride<4> > getThetas() const
    {
        return Eigen::Map<const VecXx, Eigen::Unaligned, Eigen::InnerStride<4> >(
                m_value.data() + 3, m_numEdges );
    }

    Scalar getTheta( IndexType vtx ) const
    {
        assert( vtx < m_numEdges );

        return get()[4 * vtx + 3];
    }

    void setThetas( const VecXx& thetas, int numberOfFixedThetas = 0 )
    {
        assert( thetas.size()==m_numEdges );

        Eigen::Map<VecXx, Eigen::Unaligned, Eigen::InnerStride<4> >(
                m_value.data() + 4 * numberOfFixedThetas + 3, m_numEdges - numberOfFixedThetas ) =
                thetas.tail( m_numEdges - numberOfFixedThetas );
        setDependentsDirty();
    }

    void setTheta( IndexType vtx, Scalar theta )
    {
        m_value[4 * vtx + 3] = theta;
        setDependentsDirty();
    }

    IndexType getNumEdges() const
    {
        return m_numEdges;
    }

    IndexType getNumVertices() const
    {
        return m_numEdges + 1;
    }

    virtual const char* name() const
    {
        return "DOFs";
    }

protected:
    virtual void compute() // Not implemented as this is an pure input node
    {
        ErrorStream( g_log, "" ) << "DegreesOfFreedom::compute() should never be called";
    }

private:
    IndexType m_numEdges;
};

/**
 * Unit: cm
 */
class Edges: public DependencyNode<Vec3xArray>
{
public:
    Edges( DOFs& dofs ) :
            DependencyNode<Vec3xArray>( 0, dofs.getNumEdges() ), 
            m_dofs( dofs )
    {
        m_dofs.addDependent( this );
    }

    virtual const char* name() const
    {
        return "Edges";
    }

protected:
    virtual void compute();

    DOFs& m_dofs;
};

/**
 * Unit: cm
 */
class Lengths: public DependencyNode<std::vector<Scalar> >
{
public:
    Lengths( Edges& edges ) :
            DependencyNode<std::vector<Scalar> >( 0, edges.size() ), 
            m_edges( edges )
    {
        m_edges.addDependent( this );
    }

    virtual const char* name() const
    {
        return "Lengths";
    }

protected:
    virtual void compute();

    Edges& m_edges;
};

/**
 * Unit: no dimension
 */
class Tangents: public DependencyNode<Vec3xArray>
{
public:
    Tangents( Edges& edges, Lengths& lengths ) :
            DependencyNode<Vec3xArray>( 0, edges.size() ), 
            m_edges( edges ), 
            m_lengths( lengths )
    {
        m_edges.addDependent( this );
        m_lengths.addDependent( this );
    }

    virtual const char* name() const
    {
        return "Tangents";
    }

protected:
    virtual void compute();

    Edges& m_edges;
    Lengths& m_lengths;
};

/**
 * Unit: no dimension
 */
class CurvatureBinormals: public DependencyNode<Vec3xArray>
{
public:
    CurvatureBinormals( Tangents& tangents ) :
            DependencyNode<Vec3xArray>( 1, tangents.size() ), 
            m_tangents( tangents )
    {
        m_tangents.addDependent( this );
    }

    virtual const char* name() const
    {
        return "CurvatureBinormals";
    }

protected:
    virtual void compute();

    Tangents& m_tangents;
};

/**
 * Unit: no dimension
 */
class TrigThetas: public DependencyNode<std::pair<VecXx, VecXx> >
{
public:
    TrigThetas( DOFs& dofs ) :
            DependencyNode<std::pair<VecXx, VecXx> >( std::make_pair( VecXx(), VecXx() ) ), 
            m_dofs( dofs )
    {
        m_dofs.addDependent( this );
    }

    virtual const char* name() const
    {
        return "TrigThetas";
    }

    const VecXx& getSines()
    {
        return get().first;
    }

    const VecXx& getCosines()
    {
        return get().second;
    }

protected:
    virtual void compute();

    DOFs& m_dofs;

#ifndef WETA
private:
    void vdSinCos( const int n, const double a[], double r1[], double r2[] );
#endif

};

#endif /* DEGREESOFFREEDOM_HH_ */
