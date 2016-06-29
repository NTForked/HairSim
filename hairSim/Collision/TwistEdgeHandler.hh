#ifndef TWIST_EDGE_HANDLER_H
#define TWIST_EDGE_HANDLER_H

#include "ElementProxy.hh"
#include "../Core/ElasticStrand.hh"
#include "../Dynamic/ImplicitStepper.hh"

class TwistEdge;

namespace strandsim
{
    class TwistEdgeHandler 
    {
    public:

        TwistEdgeHandler();
        ~TwistEdgeHandler();

        void getEdgeVerts( const TwistEdge* edge, const bool& startOfStep, Vec3x& firstVert, Vec3x& secondVert );
        void getEdgeVerts( const TwistEdge* edge, const bool& startOfStep, Vec3x& firstVert, Vec3x& secondVert, const VecXx& vels, const double& dt );
        void getEdgeAlphas( const TwistEdge* edge, const bool& startOfStep, double& alpha, double& beta );

        void updateTwistAngle( TwistEdge* edge, TwistEdge* startPA, TwistEdge* startPB, TwistEdge* endPA, TwistEdge* endPB, const bool& traversal );
        void deleteBand( TwistEdge* edge, std::vector< ElementProxy* >& elementProxies );
        void applyImpulses( ElasticStrand* strand, ImplicitStepper* stepper, bool computeImpulses );

        std::vector< TwistEdge* > tunnelingBands;

        int m_frozenCheck;
        bool m_frozenScene;
        bool trackTunneling;
        bool repeatCD;
    };
}
#endif

