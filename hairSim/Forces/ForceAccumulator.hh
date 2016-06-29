#ifndef FORCEACCUMULATOR_HH_
#define FORCEACCUMULATOR_HH_

#include "../Utils/Definitions.h"
#include "../Math/BandMatrix.h"

class ElasticStrand;
class StrandState;

template<typename ForceT> class ForceAccumulator
{
public:
    template<typename AccumulatedT>
    static void accumulateCurrent( AccumulatedT& accumulated, ElasticStrand& strand )
    {
        // Yes, this statement does nothing... It is there to guarantee that accumulateCurrent
        // will never be called on dissipative forces, as they require evolution to be estimated.
        // So if you get a compilation error here, that's why. If thas is a problem just comment out
        // the line below, but make sure that dissipative forces still make sense when called on
        // the current state (most likely they should return zero).
        ForceT::NonDissipativeForce;

        accumulate( accumulated, strand, strand.getCurrentState() );
    }

    template<typename AccumulatedT>
    static void accumulateFuture( AccumulatedT& accumulated, ElasticStrand& strand )
    {
        accumulate( accumulated, strand, strand.getFutureState() );
    }

    static void accumulate( Scalar& energy, const ElasticStrand& strand, StrandState& state )
    {
        for ( IndexType vtx = ForceT::s_first; vtx < state.m_numVertices - ForceT::s_last; ++vtx )
        {
            energy += ForceT::localEnergy( strand, state, vtx );
        }
    }

    static void accumulate( VecXx& force, const ElasticStrand& strand, StrandState& state )
    {
        typename ForceT::LocalForceType localF;
        for ( IndexType vtx = ForceT::s_first; vtx < state.m_numVertices - ForceT::s_last; ++vtx )
        {
            ForceT::computeLocal( localF, strand, state, vtx );
            ForceT::addInPosition( force, vtx, localF );
        }
    }

    static void accumulate( JacobianMatrixType& Jacobian, const ElasticStrand& strand, StrandState& state )
    {
        typename ForceT::LocalJacobianType localJ;
        for ( IndexType vtx = ForceT::s_first; vtx < state.m_numVertices - ForceT::s_last; ++vtx )
        {
            ForceT::computeLocal( localJ, strand, state, vtx );
            ForceT::addInPosition( Jacobian, vtx, localJ );
        }
    }

};

#endif /* FORCEACCUMULATOR_HH_ */
