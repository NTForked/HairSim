#ifndef SIMULATION_H
#define SIMULATION_H

#include "SimulationParameters.hh"
#include "ProblemFwd.hh"

#include "../Core/Definitions.hh"
#include "../Utils/SpatialHashMapFwd.hh"
#include "../Collision/ProximityCollision.hh"
#include "../Collision/EdgeEdgeCollision.hh"
#include "../../bogus/Interfaces/MecheEigenInterface.hpp"

#include <vector>
#include <list>
#include <set>
#include <map>
#include <tr1/memory>

class ElasticStrand;
class ImplicitStepper;
class CollisionDetector;
class ElementProxy;
class CollisionBase;

class Simulation
{
public:

    friend class SimulationUtils;

    Simulation( const std::vector<ElasticStrand*>& strands, Scalar& dt, const SimulationParameters& params );

    virtual ~Simulation();

    //! Step the simulation forward by m_dt
    /*! \param substepsInDt the number of simulation substeps to perform */
    void execute( int substepsInDt );

    void setDt( Scalar dt )
    { m_dt = dt; }

    Scalar getDt() const
    { return m_dt; }

private:

    //! Performs a simulation substep
    void step( unsigned subStepId, Scalar dt );

    bool isCollisionInvariantCT( const Scalar dt );
    void deleteInvertedProxies();
    void traversalCheck();

    //! Prepares the substep, executes the externals objects controllers
    void step_prepare( Scalar dt );
    //! Solves the unconstrained dynamics of each strand
    void step_dynamics( unsigned subStepId, Scalar dt );
    //! Analyses the current collisions and compute the colliding groups
    void step_processCollisions( Scalar dt );
    //! Solves each colliding group
    void step_solveCollisions();


//// put below funcs into separate CollisionHandler Class

    //! Returns wether a strand needs to be solved using YacFS/bogus.
    /*! Will be true if the strand is subject to at least one contact or hard constraint */
    bool needsExternalSolve( unsigned strandIdx ) const;
                
    bool assembleBogusFrictionProblem( CollidingGroup& collisionGroup, bogus::MecheFrictionProblem& mecheProblem,
            std::vector<unsigned> &globalIds, std::vector<ProximityCollision*> &colPointers, VecXx& vels, VecXx& impulses,
            VecXu& startDofs, VecXu& nDofs, int& numSubSys );
    //! Cleanup a friction problem and updates the strands with the new velocities if \p accept is true
    void postProcessBogusFrictionProblem( bool accept, CollidingGroup& collisionGroup,
            const bogus::MecheFrictionProblem& mecheProblem, const std::vector<unsigned> &globalIds,
            const std::vector<ProximityCollision*> &colPointers, VecXx& vels, VecXx& impulses,
            VecXu& startDofs, VecXu& nDofs  );
    //! Proper solving of the MecheFrictionProblem
    bool solveBogusFrictionProblem( bogus::MecheFrictionProblem& mecheProblem, const std::vector<unsigned> &globalIds,
            bool asFailSafe, bool nonLinear, VecXx& vels, VecXx& impulses, int& numSubSys );

    //! Solve the contacts and constraints on a single object
    void solveSingleObject( const unsigned objectIdx, bool asFailSafe, bool nonLinear );
    //! Solve the contacts and constraints on a colliding group
    void solveCollidingGroup( CollidingGroup &cg, bool asFailSafe, bool nonLinear );

    //! previously Mesh/hair collision detection
    void setupContinuousTimeCollisions();
    //! previously only Hair/hair collision detection
    void setupProximityCollisions( Scalar dt );

    //! Continuous-time mesh/hair collisision detection
    void doContinuousTimeDetection( Scalar dt );

    //! Adds an external contact on strand \p strIdx, edge \p edgeIdx, abscissa \p abscissa
    /*! \return whether this collision has been accepted */
    bool addExternalContact( const unsigned strIdx, const unsigned edgeIdx, const Scalar abscissa, const ProximityCollision& collision );

    //! Computes the deformation gradient of a strand at one contact point, ie dq/dx
    void computeDeformationGradient( ProximityCollision::Object &object ) const;
    //! Setup the local frame for one contact and calls computeDeformationGradient() for each object
    void setupDeformationBasis( ProximityCollision &collision ) const;

    //! Transform a mutual collision into an external contact on the (onFirstObject ? first : second) object
    void makeExternalContact( ProximityCollision& c, bool onFirstObject );

    //! Computes the colliding groups using a graph walking algorithm
    void computeCollidingGroups( const ProximityCollisions &mutualCollisions, Scalar dt );

    //! Returns whether a collision is deemed acceptable ( not too close to the root, etc )
    static bool acceptsCollision( const ElasticStrand& strand, int edgeIdx, Scalar localAbscissa,
            const Vec3x& normal );


//////////

    Scalar& m_dt;   //!< Time per "frame"
    SimulationParameters& m_params; // there should only be one, this is a reference to Scene's simParams
    const std::vector< ElasticStrand* >& m_strands;

    std::vector< ImplicitStepper* > m_steppers;
    CollisionDetector* m_collisionDetector;        //!< BVH-based collision detector

    //! Map between a index in the simulation to an index in a colliding group
    typedef std::map<unsigned, unsigned> IndicesMap;
    //! Colliding group: set of strands and contacts that should be solved together
    typedef std::pair<IndicesMap, ProximityCollisions> CollidingGroup;


    std::vector<ProximityCollisions> m_externalContacts;  //!< External contacts on each strand
    ProximityCollisions m_mutualContacts;           //!< List of all rod-rod contacts

    std::vector<CollidingGroup> m_collidingGroups;
    std::vector<unsigned> m_globalIds;

    //! Index of colliding group in which each strand should be. Can be -1.
    std::vector<int> m_collidingGroupsIdx;

    //!< Spatial Hash Map for hair/hair proximity collision detetection
    typedef SpatialHashMap<ElasticStrand, unsigned, true> SpatialHashMapT;
    SpatialHashMapT * m_hashMap;

    int m_num_ct_hair_hair_col;
};

#endif
