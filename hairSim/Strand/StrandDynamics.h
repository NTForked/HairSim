#ifndef STRANDSIM_STRANDDYNAMICTRAITS_HH
#define STRANDSIM_STRANDDYNAMICTRAITS_HH

#include "../Core/Definitions.hh"
#include "../Core/BandMatrixFwd.hh"

class ElasticStrand ;
class DOFScriptingController;

class StrandDynamicTraits
{

public:
    StrandDynamicTraits( ElasticStrand& strand );
    ~StrandDynamicTraits() ;

    void resizeSelf() ;

    VecXx& getDisplacements()
    {
        this needs to get reworked where used, in order to be able to modify these displacements
        return m_strand.getFutureDegreesOfFreedom() - m_strand.getCurrentDegreesOfFreedom()
    }

    void setDisplacements( const Vec3x& difference )
    {
        m_strand.setFutureDegreesOfFreedom( m_strand.getCurrentDegreesOfFreedom() + difference );
    }

    Vec3x getDisplacement( IndexType vtx ) const
    {
        return m_strand.getFutureVertex( vtx ) - m_strand.getVertex( vtx );
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
    void acceptGuess();

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
    DOFScriptingController* m_scriptingController;

};

#endif // STRANDSIM_STRANDDYNAMICTRAITS_HH
