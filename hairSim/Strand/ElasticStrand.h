#ifndef ELASTICSTRAND_HH_
#define ELASTICSTRAND_HH_

#include <list>
#include <set>
#include <tr1/memory>

#include "ElasticStrandParameters.h"

#include "../Forces/ForceBase.hh"

#include "StrandBase.h"
#include "StrandState.h"
#include "BandMatrixFwd.h"
#include "CollisionParameters.hh"
#include "../Utils/ThreadUtils.hh"
#include "../../Apps/StrandSimulator/ProblemStepper.hh"

template<typename ViscousT> class StretchingForce;
template<typename ViscousT> class BendingForce;
template<typename ViscousT> class TwistingForce;
template<typename ForceT> class ForceAccumulator;

class GravitationForce;
class AirDragForce;
class CollisionSet;
class StrandDynamicTraits;
class DOFScriptingController;

class ElasticStrand
{
public:

    ElasticStrand( const VecXx& dofs, const ElasticStrandParameters& parameters,
            DOFScriptingController* controller, double collisionRadius = 0.01, 
            int globalIndex = -1, const Vec3x& initRefFrame1 = Vec3x() );

    virtual ~ElasticStrand();

    IndexType getNumVertices() const { return m_numVertices; }
    IndexType getNumEdges() const { return m_numEdges; }

    // Acknowledges that current physics are not valid anymore. Which happens more often than one may think.
    // This could be due to several things, for instance a change in physical parameters or external forces.
    // For statics, it means that we may no longer be in an equilibrium positions ; we have to resume updating
    void invalidatePhysics();
    // Invalidate physics of future state only
    void invalidateFuturePhysics();

    void swapStates()
    {
        // std::cout << "SWAP [" << m_globalIndex << "]" <<  std::endl;
        std::swap( m_currentState, m_futureState );
    }

    const JacobianMatrixType& getTotalJacobian() const
    {
        return m_totalJacobian;
    }

    JacobianMatrixType& getTotalJacobian()
    {
        return m_totalJacobian;
    }


    ////////////////////////////////////////////////////////////////////////////////
    // Access to current state

    const VecXx& getCurrentDegreesOfFreedom() const
    {
        return m_currentState->m_dofs.get();
    }
    void setCurrentDegreesOfFreedom( const VecXx& dof );

    Vec3x getVertex( IndexType vtx ) const
    {
        return m_currentState->getVertex( vtx );
    }
    void setVertex( IndexType vtx, const Vec3x& point )
    {
        m_currentState->setVertex( vtx, point );
    }
    Scalar getTheta( const IndexType vtx ) const
    {
        assert( vtx < m_numEdges );
        return m_currentState->getTheta( vtx );
    }

    void setTheta( const IndexType vtx, const Scalar newTheta )
    {
        m_currentState->setTheta( vtx, newTheta );
        invalidatePhysics();
    }

    const Vec3x& getCurrentTangent( int vtx ) const
    {
        return m_currentState->m_tangents[vtx];
    }
    const Vec3xArray& getCurrentReferenceFrames1()
    {
        return m_currentState->m_referenceFrames1.get();
    }
    const Vec3xArray& getCurrentReferenceFrames2()
    {
        return m_currentState->m_referenceFrames2.get();
    }
    const Vec3xArray& getCurrentMaterialFrames1()
    {
        return m_currentState->m_materialFrames1.get();
    }
    const Vec3xArray& getCurrentMaterialFrames2()
    {
        return m_currentState->m_materialFrames2.get();
    }
    const std::vector<Scalar>& getCurrentReferenceTwists() const
    {
        return m_currentState->m_referenceTwists.get();
    }

    const std::vector<Scalar>& getCurrentReferenceTwistsDirty() const
    {
        return m_currentState->m_referenceTwists.getDirty();
    }
    const std::vector<Scalar>& getFutureReferenceTwistsDirty() const
    {
        return m_futureState->m_referenceTwists.getDirty();
    }

    void setCurrentReferenceFrames1( const Vec3xArray& reff1 )
    {
        m_currentState->m_referenceFrames1.set( reff1 );
    }
    void setCurrentReferenceFrames2( const Vec3xArray& reff2 )
    {
        m_currentState->m_referenceFrames2.set( reff2 );
    }
    void setCurrentReferenceTwists( const std::vector<Scalar>& reft )
    {
        m_currentState->m_referenceTwists.set( reft );
    }

    void setCurrentReferenceTwistsClean( const std::vector<Scalar>& reft )
    {
        m_currentState->m_referenceTwists.cleanSet( reft );
    }

    void setFutureReferenceTwistsClean( const std::vector<Scalar>& reft )
    {
        m_futureState->m_referenceTwists.cleanSet( reft );
    }

    Vec3x getMaterialFrame1( int vtx ) const
    {
        return m_currentState->getMaterialFrame1( vtx );
    }

    Vec3x getMaterialFrame2( int vtx ) const
    {
        return m_currentState->getMaterialFrame2( vtx );
    }

    Vec3x getEdgeVector( int vtx ) const
    {
        return m_currentState->getEdgeVector( vtx );
    }

    Scalar getTotalEnergy() const
    {
        return m_currentState->m_totalEnergy;
    }

    const VecXx& getTotalForces() const
    {
        return m_currentState->m_totalForce;
    }

    VecXx& getTotalForces()
    {
        return m_currentState->m_totalForce;
    }

    Scalar getUnsignedAngleToMajorRadius( int vtx, const Vec3x& vec ) const;
    Scalar getSignedAngleToMajorRadius( int vtx, const Vec3x& vec ) const;
    Scalar getCurrentTotalLength() const;
    void freezeRestShape( unsigned begin, unsigned end, Scalar damping = 0. );

    ///////////////////////////////////////////
    // Future geometry

    Vec3x getFutureVertex( IndexType vtx ) const
    {
        return m_futureState->getVertex( vtx );
    }

    const VecXx& getFutureDegreesOfFreedom() const
    {
        return m_futureState->m_dofs.get();
    }
    void setFutureDegreesOfFreedom( const VecXx& dof );

    Scalar getFutureTotalLength() const;

    void filterFutureGeometryByRestLength( const double epsilon = 0., bool allowCompression = false );

    VecXx& getFutureTotalForces()
    {
        return m_futureState->m_totalForce;
    }

    Scalar getFutureTotalEnergy() const
    {
        return m_futureState->m_totalEnergy;
    }

    const VecXx& getFutureTotalForces() const
    {
        return m_futureState->m_totalForce;
    }

    VecXx& getFutureTotalForces()
    {
        return m_futureState->m_totalForce;
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Changing the shape

    void findBendElasticLimitMax( Scalar& bendElasticLimit ) const;

    void applyPlasticDeformation( Scalar stretchElasticLimit, Scalar bendElasticLimit,
            Scalar twistElasticLimit );

    void updateEverythingThatDependsOnRestLengths();

    /////////////////////////////
    // Rest Shape and physical parameters

    Scalar getEdgeRestLength( const IndexType vtx ) const
    {
        assert( vtx < m_numEdges );
        return m_restLengths[vtx];
    }

    void setEdgeRestLength( const IndexType vtx, const Scalar newrestlength );

    void setEdgesRestLength( const Scalar newRestLength );

    void setRestShape( const VecXx &dofs, unsigned begin, unsigned end, Scalar damping = 0. );

    Scalar getRadiusA( int vtx ) const
    {
        return m_parameters.getRadiusA( vtx );
    }

    Scalar getRadiusB( int vtx ) const
    {
        return m_parameters.getRadiusB( vtx );
    }

    void setRadius( const Scalar radius_a, const Scalar radius_b );

    Scalar getRestTwist( const IndexType vtx ) const
    {
        return m_restTwists[vtx];
    }

    void setRestTwist( const IndexType vtx, const Scalar restTwist )
    {
        m_restTwists[vtx] = restTwist;
        invalidatePhysics();
    }

    Vec2x getKappaBar( const IndexType vtx ) const
    {
        return m_restKappas[vtx];
    }

    void setKappaBar( const IndexType vtx, const Vec2x& kappaBar )
    {
        m_restKappas[vtx] = kappaBar;
        invalidatePhysics();
    }

    Scalar getStiffness() const
    {
        return m_parameters.getYoungsModulus();
    }

    void setStiffness( const Scalar youngs );

    void setParameters( double i_radiusA, double i_radiusB, double i_rootRM, double i_tipRM,
            double i_youngsModulus, double i_shearModulus, double i_density, double i_viscosity,
            double i_airDrag );

    void setParameters( const ElasticStrandParameters &parameters );

    const ElasticStrandParameters &getParameters() const
    {
        return m_parameters;
    }
    ElasticStrandParameters &getParameters()
    {
        return m_parameters;
    }

    //! Informs the strand that its parameters have been modified
    void ackParametersChanged();

    CollisionParameters& collisionParameters()
    {
        return m_collisionParameters;
    }

    const CollisionParameters& collisionParameters() const
    {
        return m_collisionParameters;
    }

    Scalar getTotalRestLength() const
    {
        return m_totalRestLength;
    }

    Scalar getVertexMass( IndexType i ) const
    {
        return m_vertexMasses[i];
    }

    Scalar getEdgeInertia( IndexType i )
    {
        const Scalar a = m_parameters.getRadiusA( i );
        const Scalar b = m_parameters.getRadiusB( i );
        const Scalar mass = m_parameters.getDensity() * M_PI * a * b * m_restLengths[i];

        return 0.25 * mass * ( square( a ) + square( b ) );
    }

    Scalar getCurvilinearAbscissa( int vtx, Scalar localAbscissa ) const;
    void getLocalAbscissa( const Scalar curvilinearAbscissa, int &vtx,
            Scalar &localAbscissa ) const;

    //////////////////////////////////////////////
    // States, dynamics

    StrandState& getCurrentState()
    {
        return *m_currentState;
    }

    const StrandState& getCurrentState() const
    {
        return *m_currentState;
    }

    const StrandState& getFutureState() const
    {
        return *m_futureState;
    }

    StrandState& getFutureState()
    {
        return *m_futureState;
    }

    bool canDoDynamics() const
    {
        return m_dynamics != NULL;
    }
    void createDynamics();

    const StrandDynamics& dynamics() const
    {
        return *m_dynamics;
    }

    StrandDynamics& dynamics()
    {
        return *m_dynamics;
    }

    bool requiresExactJacobian() const
    {
        return m_requiresExactJacobian;
    }

    void requireExactJacobian( bool b )
    {
        m_requiresExactJacobian = b;
    }

    bool activelySimulated() const
    {
        return m_activelySimulated;
    }

    void setActivelySimulated( bool simulated )
    {
        m_activelySimulated = simulated;
    }

    // Operates on future state's thetas
    void relaxThetas();
    void resetThetas();

    ///////////////
    // Segments 

    void getSegment( unsigned elementID, Vec3x &start, Vec3x &end ) const;

    // Begin and end of edge iterators, used in SpatialHashMap
    unsigned subsamples_begin() const
    { return 0; }

    unsigned subsamples_end() const
    { return m_numEdges; }

    //////////////////////
    // Serialization

    bool serializeTo( std::ostream & os ) const;
    bool deserializeFrom( std::istream & is );

    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
    {
        ar & m_currentState->m_dofs;
        ar & m_currentState->m_referenceFrames1;
        ar & m_currentState->m_referenceFrames2;
        ar & m_currentState->m_referenceTwists;
        ar & m_currentState->m_referenceFrames1.getPreviousTangents();
        ar & m_parameters.dt();
    }

    MutexType& mutex()
    {
        return *m_mutex;
    }

    int getGlobalIndex() const
    {
        return m_globalIndex;
    }

    void setGlobalIndex( int globalIndex )
    {
        m_globalIndex = globalIndex;
    }

    const Vec2xArray& getRestKappas() const
    {
        return m_restKappas;
    }

    Vec2xArray& alterRestKappas()
    {
        return m_restKappas;
    }

    const std::vector<Scalar>& getRestLengths() const
    {
        return m_restLengths;
    }

    const std::vector<Scalar>& getRestTwists() const
    {
        return m_restTwists;
    }

private:

    void resizeInternals();

    /////////////////////////////////////
    // Convenience force accumulators

    // Add ForceT's energy to the geometry
    template<typename ForceT>
    void accumulateE( StrandState& geometry ) const
    {
        ForceAccumulator<ForceT>::accumulate( geometry.m_totalEnergy, *this, geometry );
    }

    // Add ForceT to the geometry
    template<typename ForceT>
    void accumulateF( StrandState& geometry ) const
    {
        ForceAccumulator<ForceT>::accumulate( geometry.m_totalForce, *this, geometry );
    }

    // Add ForceT's Jacobian to the geometry
    template<typename ForceT>
    void accumulateJ( StrandState& geometry ) const
    {
        ForceAccumulator<ForceT>::accumulate( *geometry.m_totalJacobian, *this, geometry );
    }

    // Add ForceT's energy and force to the geometry
    template<typename ForceT>
    void accumulateEF( StrandState& geometry ) const
    {
        accumulateE<ForceT>( geometry );
        accumulateF<ForceT>( geometry );
    }

    // Add ForceT's energy, force and Jacobian to the geometry
    template<typename ForceT>
    void accumulateEFJ( StrandState& geometry ) const
    {
        accumulateE<ForceT>( geometry );
        accumulateF<ForceT>( geometry );
        accumulateJ<ForceT>( geometry );
    }
    // Add ForceT's energy to the geometry
    template<typename ForceT>
    void accumulate( ForceBase::Quantities q, StrandState& geometry ) const
    {
        switch ( q )
        {
        case ForceBase::EFJ:
            accumulateJ<ForceT>( geometry );
            /* no break */
        case ForceBase::EF:
            accumulateF<ForceT>( geometry );
            /* no break */
        case ForceBase::E:
            accumulateE<ForceT>( geometry );
            break;
        case ForceBase::F:
            accumulateF<ForceT>( geometry );
            break;
        case ForceBase::J:
            accumulateJ<ForceT>( geometry );
            break;
        case ForceBase::NONE:
            break;
        };
    }

    template<typename ForceT>
    void accumulateEFThetaOnly( Scalar& thetaEnergy, VecXx& thetaForce,
            StrandState& geometry ) const;
    template<typename ForceT>
    void accumulateJThetaOnly( TriDiagonalMatrixType& J, StrandState& geometry ) const;

    //////////////////////////////////////////////

    double m_collisionRadius;
    double m_physicsRadius;

    int m_globalIndex; // Global index in the simulation

    // Size of the strand. The original belongs to m_parameters, so it can correctly interpolate when asked e.g. for a variable radius.
    IndexType& m_numVertices;
    IndexType m_numEdges;

    // Other physical parameters
    ElasticStrandParameters m_parameters;
    CollisionParameters m_collisionParameters;

    // Current and future geometry. They need to be pointers to be easily swapped.
    StrandState* m_currentState;
    StrandState* m_futureState;

    // Dynamics & controller
    StrandDynamics* m_dynamics;

    // The Jacobian is only needed once. To avoid memory duplication, share between current and future geometries.
    JacobianMatrixType m_totalJacobian;

    // Rest shape
    std::vector<Scalar> m_restLengths; // The following four members depend on m_restLengths, which is why updateEverythingThatDependsOnRestLengths() must be called
    Scalar m_totalRestLength;
    std::vector<Scalar> m_VoronoiLengths; // rest length around each vertex
    std::vector<Scalar> m_invVoronoiLengths; // their inverses
    std::vector<Scalar> m_vertexMasses;
    Vec2xArray m_restKappas;
    std::vector<Scalar> m_restTwists;

    // Flags
    bool m_requiresExactJacobian;
    bool m_activelySimulated;

    // Other stuff
    MutexWrapper m_mutex;

    friend class Viscous;
    friend class NonViscous;
    template<typename ViscousT> friend class StretchingForce;
    template<typename ViscousT> friend class BendingForce;
    template<typename ViscousT> friend class TwistingForce;
    friend class GravitationForce;
    friend class AirDragForce;
    friend class StrandDynamics;
    friend std::ostream& operator<<( std::ostream& os, const ElasticStrand& strand );

};


BOOST_CLASS_VERSION( strandsim::ElasticStrand, 0 )

#endif /* ELASTICSTRAND_HH_ */
