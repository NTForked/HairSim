#ifndef TWIST_EDGE_HANDLER_H
#define TWIST_EDGE_HANDLER_H

#include "ElementProxy.hh"
#include "../Core/ElasticStrand.hh"
#include "../Dynamic/ImplicitStepper.hh"

namespace strandsim
{

    struct TwistEdge;
    
    class TwistEdgeHandler 
    {
    public:

        TwistEdgeHandler( const double& thickness);
        ~TwistEdgeHandler();

        void getEdgeVerts( const TwistEdge* edge, const bool& startOfStep, Vec3x& firstVert, Vec3x& secondVert );
        void getEdgeVerts( const TwistEdge* edge, const bool& startOfStep, Vec3x& firstVert, Vec3x& secondVert, const VecXx& vels, const double& dt );
        void getEdgeAlphas( const TwistEdge* edge, const bool& startOfStep, double& alpha, double& beta );
        // bool assignTunnSign( TwistEdge* first, TwistEdge* second ); // dont need
        bool computeTunnSign( TwistEdge* first, TwistEdge* second, const bool& startOfStep );
        // bool computeTunnSign( TwistEdge* edge ); // dont need
        // void updateTwistAngle( TwistEdge* edge );
        void updateTwistAngle( TwistEdge* edge, TwistEdge* startPA, TwistEdge* startPB, TwistEdge* endPA, TwistEdge* endPB, const bool& traversal );
        void deleteBand( TwistEdge* edge, std::vector< ElementProxy* >& elementProxies );
        void applyImpulses( ElasticStrand* strand, ImplicitStepper* stepper, bool computeImpulses );

        std::vector< TwistEdge* > tunnelingBands;
        std::pair< int, std::pair< Vec3x, Vec3x > > qs_prime;
        std::pair< int, std::pair< Vec3x, Vec3x > > qe_prime;

        int m_frozenCheck;
        bool m_frozenScene;
        bool trackTunneling;
        bool repeatCD;
        const double m_thickness;
    };
}
#endif

