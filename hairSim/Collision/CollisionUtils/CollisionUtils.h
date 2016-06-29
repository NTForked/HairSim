#ifndef COLLISION_UTILS_H
#define COLLISION_UTILS_H

#include "../../Utils/Definitions.h"
#include "../CollisionParameters.h"
#include "../Collision.h"

class ElasticStrand;
class StrandState;

/////////////////////////////////////////////////////////////////
// Some of this is Collision detection code adapted from Robert Bridson's website

void addUnique( std::vector<double>& a, double e );

double triple( const Vec3 &a, const Vec3 &b, const Vec3 &c );

double signed_volume( const Vec3& x0, const Vec3& x1, const Vec3& x2, const Vec3& x3 );

void getCoplanarityTimes( 
		const Vec3& x0, const Vec3& x1, const Vec3& x2, const Vec3& x3,
        const Vec3& xnew0, const Vec3& xnew1, const Vec3& xnew2, const Vec3& xnew3,
        double *times, double *errors, unsigned &num_times );

void getIntersectionPoint( const Vec3& x0, const Vec3& xnew0, const Vec3& xnew1,
        const Vec3& xnew2, const Vec3& xnew3, double *times, double *errors,
        unsigned &num_times );

bool analyseRoughRodRodCollision( const ElasticStrand* sP, const ElasticStrand* sQ, const int iP,
        const int iQ, Vec3 &normalQtoP, Scalar &s, Scalar &t, Scalar &distance );

bool compareCT( const Collision* ct1, const Collision* ct2 );

#endif

// TODO:
//     It would be nice to handle degenerate cases better in these methods.  
//     They all handle degenerate cases, but getting PREDICTABLE behavior out would rock!!!
