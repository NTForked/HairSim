#ifndef SIMULATIONPARAMETERS_HH_
#define SIMULATIONPARAMETERS_HH_

#include "../Utils/Definitions.h"

struct SimulationParameters
{
    SimulationParameters():
        m_useProxRodRodCollisions( true ),
        m_useCTRodRodCollisions( false ),
        m_alwaysUseNonLinear( true ),
        m_useLengthProjection( false ),
        m_inextensibility_threshold( 1. ),
        m_stretching_threshold( 2.0 ),
        m_costretch_residual_threshold( 0.0 ),
        m_stretchDamping( 0. )
    {}

    int m_numberOfThreads;
    bool m_simulationManager_limitedMemory;
    unsigned m_maxNewtonIterations;

    /**
     * Rod-rod collisions
     */
    
    bool m_useProxRodRodCollisions; // whether we should use rod-rod proximity collisions (TODO: should rename; keeping name for now to maintain backwards comaptibility of older examples)
    bool m_useCTRodRodCollisions; // whether we should use rod-rod ctc collisions

    bool m_useNonLinearAsFailsafe;
    bool m_alwaysUseNonLinear;
    
    /**
     * Inextensibility 
     */
    bool m_useLengthProjection;
    double m_inextensibility_threshold; // (not in use) # of times the original step has to be halved before the inextensibility filter is applied
    double m_stretching_threshold;
    double m_costretch_residual_threshold;
    double m_stretchDamping;
    double m_gaussSeidelTolerance;

};

#endif 
