#ifndef STRAND_DYNAMICTRAITS_H
#define STRAND_DYNAMICTRAITS_H

#include "../Utils/Definitions.h"
#include "../Math/BandMatrixFwd.h"
#include "ElasticStrand.h"

class DOFScriptingController;

class StrandDynamics
{

public:
    StrandDynamics( ElasticStrand& strand );
    ~StrandDynamics() ;

    void resizeSelf() ;

    const VecXx getDisplacements() const
    {
        return m_strand.getFutureDegreesOfFreedom() - m_strand.getCurrentDegreesOfFreedom();
    }

    void setDisplacements( const VecXx& difference )
    {
        m_strand.setFutureDegreesOfFreedom( m_strand.getCurrentDegreesOfFreedom() + difference );
    }

    Vec3 getDisplacement( IndexType vtx ) const
    {
        return m_strand.getFutureVertex( vtx ) - m_strand.getVertex( vtx );
    }

    VecXx getFutureVelocities( const Scalar dt ) const
    {
        return getDisplacements() / dt;
    }

    VecXx getCurrentVelocities() const
    {
        return m_currentVelocities; // start of step vel (finite difference previous step displacements / previous step dt)
    }

    void invalidatePhysics()
    {
        m_futureForcesUpToDate = false ;
        m_futureJacobianUpToDate = false ;
        m_DOFmassesUpToDate = false ;
    }
    void invalidateFuturePhysics()
    {
        m_futureForcesUpToDate = false ;
        m_futureJacobianUpToDate = false ;
    }

    // Flags
    bool isFutureJacobianUpToDate() const
    { return m_futureJacobianUpToDate; }

    void setFutureJacobianUpToDate( bool ok )
    { m_futureJacobianUpToDate = ok; }

    bool isFutureForcesUpToDate() const
    { return m_futureForcesUpToDate; }


    // Dynamic
    const VecXx& getDOFMasses() const;
    void computeDOFMasses();
    void computeViscousForceCoefficients(Scalar dt) ;
    void computeFutureJacobian( bool withViscous = true, bool butOnlyForBendingModes = false );
    void computeLHS( Scalar dt, bool withViscous );
    void computeFutureForces( bool withViscous = true, bool butOnlyForBendingModes = false );
    void computeFutureConservativeEnergy();
    void addMassMatrixTo( JacobianMatrixType& J ) const;
    void multiplyByMassMatrix( VecXx& F ) const;
    void acceptFuture();

    // Controller
    void setScriptingController( DOFScriptingController *controller )
    {
        m_scriptingController = controller;
    }

    DOFScriptingController* getScriptingController()
    {
        return m_scriptingController;
    }

    /*
     * Checks if current state has NaNs; it that case replace the step with rigid motion.
     * This is supposed to be called after a step has left the pre-step position in m_futureState.
     */
    void nanFailSafe();

private:

    ElasticStrand& m_strand;

    // Flags
    bool m_futureJacobianUpToDate;
    bool m_futureForcesUpToDate;

    bool m_DOFmassesUpToDate;
    VecXx m_DOFmasses;

    VecXx m_currentVelocities; // SERIALIZE_ME for checkpointing

    DOFScriptingController* m_scriptingController;

};

#endif 
