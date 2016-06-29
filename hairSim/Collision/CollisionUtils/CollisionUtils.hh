#ifndef STRANDSIM_COLLISIONUTILS_HH
#define STRANDSIM_COLLISIONUTILS_HH

#include "../Core/Definitions.hh"
#include "../Core/CollisionParameters.hh"

class ElasticStrand;
class StrandState;

/////////////////////////////////////////////////////////////////
// Some of this is Collision detection code adapted from Robert Bridson's website

void addUnique( std::vector<double>& a, double e );

double triple( const Vec3x &a, const Vec3x &b, const Vec3x &c );

double signed_volume( const Vec3x& x0, const Vec3x& x1, const Vec3x& x2, const Vec3x& x3 );

void getCoplanarityTimes( 
		const Vec3x& x0, const Vec3x& x1, const Vec3x& x2, const Vec3x& x3,
        const Vec3x& xnew0, const Vec3x& xnew1, const Vec3x& xnew2, const Vec3x& xnew3,
        double *times, double *errors, unsigned &num_times );

void getIntersectionPoint( const Vec3x& x0, const Vec3x& xnew0, const Vec3x& xnew1,
        const Vec3x& xnew2, const Vec3x& xnew3, double *times, double *errors,
        unsigned &num_times );

void buildFrame( const Vec3x& n_hat, Vec3x& t1, Vec3x& t2 );

bool analyseRoughRodRodCollision( const ElasticStrand* sP, const ElasticStrand* sQ, const int iP,
        const int iQ, Vec3x &normalQtoP, Scalar &s, Scalar &t, Scalar &distance );

bool compareCT( const CollisionBase* ct1, const CollisionBase* ct2 );

#endif

// TODO:
//     It would be nice to handle degenerate cases better in these methods.  
//     They all handle degenerate cases, but getting PREDICTABLE behavior out would rock!!!
