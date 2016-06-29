#ifndef SIMULATION_H
#define SIMULATION_H

#include "SimulationParameters.h"

#include "../Utils/Definitions.h"
#include "../Collision/CollisionUtils/SpatialHashMapFwd.hh"
#include "../Collision/Collision.h"
#include "../Collision/EdgeEdgeCollision.h"
#include "../../bogus/Interfaces/MecheEigenInterface.hpp"

#include <vector>
#include <list>
#include <set>
#include <map>
#include <tr1/memory>

class ElasticStrand;
class ImplicitStepper;
class CollisionDetector;
class CollisionBase;
class TriMesh;
class ElementProxy;

//! Map between a index in the simulation to an index in a colliding group
typedef std::map<unsigned, unsigned> IndicesMap;
//! Colliding group: set of strands and contacts that should be solved together
typedef std::pair<IndicesMap, CollidingPairs> CollidingGroup;

class Simulation
{
public:

    friend class SimulationUtils;

    Simulation( const std::vector<ElasticStrand*>& strands, const SimulationParameters& params, const std::vector< TriMesh* >& meshes );

    virtual ~Simulation();

    //! Performs a simulation substep
    void step( const Scalar& dt );

private:

    //! Prepares the substep, executes the externals objects controllers
    void step_prepare( Scalar dt );

    //! Solves the unconstrained dynamics of each strand
    void step_dynamics( Scalar dt );

    //! Analyses the current collisions and compute the colliding groups
    void step_processCollisions( Scalar dt );

    //! Solves each colliding group
    void step_solveCollisions();

    // take all the less important stuff out of StrandImplicitManager/Simulation and put it here
    void updateParameters( const SimulationParameters& params );

//// Sim Utils.cpp

    bool isCollisionInvariantCT( const Scalar dt );

    void accumulateProxies( std::vector< ElementProxy* > origProxys );

    void gatherProximityRodRodCollisions( Scalar dt );

    //! Returns whether a collision is deemed acceptable ( not too close to the root, etc )
    static bool acceptsCollision( const ElasticStrand& strand, int edgeIdx, Scalar localAbscissa );

    //! Transform a mutual collision into an external contact on the (onFirstObject ? first : second) object
    void makeExternalContact( CollidingPair& c, bool onFirstObject );

    void detectContinuousTimeCollisions();

    //! Continuous-time mesh/hair collisision detection
    void preProcessContinuousTimeCollisions( Scalar dt );

    //! Adds an external contact on strand \p strIdx, edge \p edgeIdx, abscissa \p abscissa
    /*! \return whether this collision has been accepted */
    bool addExternalContact( const unsigned strIdx, const unsigned edgeIdx, const Scalar abscissa, const CollidingPair& collision );

    //! Computes the colliding groups using a graph walking algorithm
    void computeCollidingGroups( const CollidingPairs &mutualCollisions );

    //! Setup the local frame for one contact and calls computeDeformationGradient() for each object
    void setupDeformationBasis( CollidingPair &collision ) const;

    //! Computes the deformation gradient of a strand at one contact point, ie dq/dx
    void computeDeformationGradient( CollidingPair::Object &object ) const;

    //! Returns wether a strand needs to be solved using YacFS/bogus.
    /*! Will be true if the strand is subject to at least one contact or hard constraint */
    bool needsExternalSolve( unsigned strandIdx ) const;    

    //! Solve the contacts and constraints on a colliding group
    void solveCollidingGroup( CollidingGroup &cg, bool asFailSafe, bool nonLinear );

    //! Solve the contacts and constraints on a single object
    void solveSingleObject( const unsigned objectIdx, bool asFailSafe, bool nonLinear );

//// SimBogusUtils.cpp

    bool assembleBogusFrictionProblem( CollidingGroup& collisionGroup, bogus::MecheFrictionProblem& mecheProblem,
            std::vector<unsigned> &globalIds, std::vector<CollidingPair*> &colPointers, VecXx& vels, VecXx& impulses,
            VecXu& startDofs, VecXu& nDofs, int& numSubSys );
    //! Cleanup a friction problem and updates the strands with the new velocities if \p accept is true
    void postProcessBogusFrictionProblem( bool accept, CollidingGroup& collisionGroup,
            const bogus::MecheFrictionProblem& mecheProblem, const std::vector<unsigned> &globalIds,
            const std::vector<CollidingPair*> &colPointers, VecXx& vels, VecXx& impulses,
            VecXu& startDofs, VecXu& nDofs  );
    //! Proper solving of the MecheFrictionProblem
    bool solveBogusFrictionProblem( bogus::MecheFrictionProblem& mecheProblem, const std::vector<unsigned> &globalIds,
            bool asFailSafe, bool nonLinear, VecXx& vels, VecXx& impulses, int& numSubSys );

//// SimTwistEdgeUtils.cpp

    void deleteInvertedProxies();

    void traversalCheck();

//// Member Variables

    SimulationParameters& m_params; // there should only be one, this is a reference to Scene's simParams
    const std::vector< ElasticStrand* >& m_strands;

    std::vector< ImplicitStepper* > m_steppers;
    CollisionDetector* m_collisionDetector; //!< BVH-based collision detector

    std::vector< CollidingPairs > m_externalContacts;  //!< External contacts on each strand
    CollidingPairs m_mutualContacts;           //!< List of all rod-rod contacts

    std::vector<CollidingGroup> m_collidingGroups;
    
    std::vector<unsigned> m_globalIds;

    //! Index of colliding group in which each strand should be. Can be -1.
    std::vector<int> m_collidingGroupsIdx;

    //!< Spatial Hash Map for hair/hair proximity collision detetection
    typedef SpatialHashMap<ElasticStrand, unsigned, true> SpatialHashMapT;
    SpatialHashMapT * m_hashMap;
};

#endif
