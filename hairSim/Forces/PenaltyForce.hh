#ifndef STRANDSIM_PENALTYFORCE_HH
#define STRANDSIM_PENALTYFORCE_HH

#include "ForceBase.hh"

#include "../Utils/IndexedSet.hh"

#include <vector>
#include <map>

class ElasticStrand;

class PenaltyForce : public TunneledBandForceBase
{
public:
    PenaltyForce();
    ~PenaltyForce();

    void setParameters( const std::vector< ElasticStrand* > &strands, const Vec3x& origin, const Scalar stiffness,
                        const Scalar restLength = 0., const bool allowCompression = false ) ;

    bool breakSprings( const std::vector< ElasticStrand* > &strands, const Scalar breakingForce ) ;

    virtual std::string getName() const
    {
        return "PenaltyForce";
    }

    virtual void accumulateEFJ( StrandState& geometry, const ElasticStrand& strand );

    virtual void accumulateEF( StrandState& geometry, const ElasticStrand& strand );

    virtual void accumulateE( StrandState& geometry, const ElasticStrand& strand );

    virtual void accumulateJ( StrandState& geometry, const ElasticStrand& strand );

    virtual void accumulateCurrentF( VecXx& globalF, const ElasticStrand& strand );

    // New style interface (consistent with static-class forces)
    virtual void accumulateFutureE( Scalar& energy, const ElasticStrand& strand );

private:

    Scalar computeLocalE( const TwistEdge* edge, const ElasticStrand& strand ) const  ;
    void computeLocalF( const TwistEdge* edge, const ElasticStrand& strand, Vec3x& localF ) const ;
    void computeLocalJ( const Vec3x &v, const ElasticStrand& strand, Mat3x &localJ ) const ;

    Scalar m_stiffness ;
    Scalar m_thickness ;
};

#endif // STRANDSIM_PenaltyForce_HH
