#ifndef DOF_SCRIPTINGCONTROLLER_HH_
#define DOF_SCRIPTINGCONTROLLER_HH_

#include "../Math/BandMatrixFwd.h"
#include <map>
#include <stack>

class StrandState;

class DOFScriptingController
{
public:
    DOFScriptingController( );
    virtual ~DOFScriptingController();

    void clear()
    {
        m_scriptedDisplacements.clear();
    }

    void freezeVertices( int vtx, bool twist = false )
    {
        int cap = 3;
        if( vtx == 0 ) twist = true;
        if( twist ) cap = 4; 
        for( int dof = 0; dof < cap; ++dof ){
            m_scriptedDisplacements[vtx * 4 + dof] = 0.0;
        }
    }
    
    void setVertexDisplacement( int vtx, const Vec3& del )
    {
        m_scriptedDisplacements[4 * vtx + 0] = del[0];
        m_scriptedDisplacements[4 * vtx + 1] = del[1];
        m_scriptedDisplacements[4 * vtx + 2] = del[2];
    }

    void setThetaDisplacement( int vtx, Scalar del )
    {
        m_scriptedDisplacements[4 * vtx + 3] = del;
    }

    void fixLHS( JacobianMatrixType& LHS ) const;
    void fixRHS( VecXx& rhs ) const;
    void fixLHSAndRHS( JacobianMatrixType& LHS, VecXx& rhs, Scalar dt ) const;

    void enforceDisplacements( VecXx& displacements ) const;
    void enforceVelocities( VecXx& velocities, Scalar dt ) const;

    /**
     * @brief Computes rigid motion based on (scripted) root vertices displacement.
     *
     * This could to be used either as initial guess in a Newton method, or as safeguard motion if
     * the regular solve fails. This DOFScriptingController generates displacements: if the first edge
     * is purely translated the translation is propagated to all DOFs; otherwise a rigid motion based
     * on parallel transport of the first edge is applied.
     *
     * @param futureDOFs: position after rigid motion
     * @param currentDOFs: position before rigid motion
     */
    void computeRigidBodyMotion( VecXx& futureDOFs, const VecXx& currentDOFs );

    /**
     * @brief Reset the strand degrees of freedom to the unsimulated positions.
     * @param rootAsWell if true, move directly the roots to the end-of-frame position and set their displacements to zero
     */
    void setToUnsimulatedPositions( VecXx& futureDOFs, const VecXx& currentDOFs, bool rootAsWell = false );

private:    
    std::map<int, Scalar> m_scriptedDisplacements;
};

#endif /* VERTEXSCRIPTINGCONTROLLER_HH_ */
