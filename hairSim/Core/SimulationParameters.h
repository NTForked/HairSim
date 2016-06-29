#ifndef SIMULATIONPARAMETERS_HH_
#define SIMULATIONPARAMETERS_HH_

#include "../Core/Definitions.h"

struct SimulationParameters
{
    SimulationParameters():
        m_rodMeshPenaltyThickness( 0.1 ), // 0.01
        m_rodMeshPenaltyStiffness( 1.0 ), //200.0
        m_rodMeshImpulsesIterations( 10 ), //
        m_rodRodImpulsesIterations( 2 ), //
        m_useProxRodRodCollisions( true ), //
        m_useCTRodRodCollisions( false ), //
        m_percentCTRodRodCollisionsAccept( 100. ), //
        m_alwaysUseNonLinear( true ),
        m_useGraphSplitOption( false ), //
        m_useLengthProjection( false ), //
        m_inextensibility_threshold( 1. ), //
        m_stretching_threshold( 2.0 ), //
        m_costretch_residual_threshold( 0.0 ), //
        m_stretchDamping( 0. ),
        m_useDeterministicSolver( true )
    {}

    int m_numberOfThreads;
    
    bool m_simulationManager_limitedMemory;
    bool m_contactProblem_limitedMemory;

    double m_rodMeshPenaltyThickness;
    double m_rodMeshPenaltyStiffness;

    int m_rodMeshImpulsesIterations;
    int m_rodRodImpulsesIterations;


    unsigned m_maxNewtonIterations ;

    /**
     * Rod-rod collisions
     */
    
    bool m_useProxRodRodCollisions; // whether we should use rod-rod proximity collisions (TODO: should rename; keeping name for now to maintain backwards comaptibility of older examples)
    bool m_useCTRodRodCollisions; // whether we should use rod-rod ctc collisions
    Scalar m_percentCTRodRodCollisionsAccept;

    bool m_useNonLinearAsFailsafe;
    bool m_alwaysUseNonLinear;

    double m_hairMeshFrictionCoefficient;
    double m_hairHairFrictionCoefficient;

    double m_stochasticPruning;
    bool m_useGraphSplitOption;
    
    /**
     * Inextensibility 
     */
    bool m_useLengthProjection;
    double m_inextensibility_threshold; // (not in use) # of times the original step has to be halved before the inextensibility filter is applied
    double m_stretching_threshold;
    double m_costretch_residual_threshold;

    double m_stretchDamping;
    bool m_useDeterministicSolver;

    double m_gaussSeidelTolerance;

};

#endif 
