#ifndef DISTANCES_HH
#define DISTANCES_HH

#include "../Utils/Definitions.h"

// Closest point on a segment to a vertex
Vec3 ClosestPtPointSegment( const Vec3& point, const Vec3& first, const Vec3& last );

Vec3 ClosestPtPointSegment( bool& extremum, const Vec3& point, const Vec3& first,
        const Vec3& last );

// Closest point on a triangle to a vertex
Vec3 ClosestPtPointTriangle( const Vec3& p, const Vec3& a, const Vec3& b, const Vec3& c );

// Computes the squared distance between and closest points of two edges.
double ClosestPtSegmentSegment( const Vec3& p1, const Vec3& q1, const Vec3& p2, const Vec3& q2,
        double& s, double& t, Vec3& c1, Vec3& c2 );

// Computes the barycentric coordiantes of a point wrt a triangle
void computeBarycentricCoordinates( const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& p,
        double& u, double& v, double& w );

// Ci: center of rectangle i ; xi and yi : extents on boths axis
// CPi: closest point on rectangle i
Scalar SquareDistRectangleToRectangle( const Vec3& C1, const Vec3 &x1, const Vec3 &y1,
        const Vec3& C2, const Vec3 &x2, const Vec3 &y2, Vec3 &CP1, Vec3 &CP2 );

// /!\ Will return false if segment and rectangle are coplanar
bool intersectionSegmentRectangle( const Vec3& s_edge_0, const Vec3& s_edge_1,
        const Vec3& r_center, const Vec3& r_extent_1, const Vec3& r_extent_2, Scalar &t );

// Templated segment/point and segment/segment distances

template<typename VecT, typename ScalarT >
ScalarT SquareDistPointToSegment( const VecT& point, const VecT& first, const VecT& last )
{
    const VecT dfirst =  ( point - first ) ;
    const ScalarT sqDistFirst = dfirst[0]*dfirst[0] + dfirst[1]*dfirst[1] + dfirst[2]*dfirst[2] ;
    if ( isSmall( sqDistFirst ) ) return sqDistFirst ;

    const ScalarT vx = last[0] - first[0];
    const ScalarT vy = last[1] - first[1];
    const ScalarT vz = last[2] - first[2];
    // Squared norm of segment
    const ScalarT len2 = ( vy * vy + vx * vx + vz * vz );
    // Abscissa of proj on line
    const ScalarT dtp = ( vx * dfirst[0] + vy * dfirst[1] + vz * dfirst[2] ) / len2 ;

    if ( dtp <= 0 )
    {
        return sqDistFirst ;
    }
    const ScalarT dtp2 = dtp*dtp ;

    if( dtp2 >= len2 )
    {
        const VecT dlast = ( point - last ) ;
        const ScalarT sqDistLast = dlast[0]*dlast[0] + dlast[1]*dlast[1] + dlast[2]*dlast[2] ;

        return sqDistLast ;
    }
    return sqDistFirst - dtp2 ;

}

template<typename VecT, typename ScalarT, typename InVecT>
ScalarT SquareDistSegmentToSegment( const Eigen::MatrixBase<InVecT>& p0,
        const Eigen::MatrixBase<InVecT>& p1, const Eigen::MatrixBase<InVecT>& q0,
        const Eigen::MatrixBase<InVecT>& q1, ScalarT &s, ScalarT &t )
{
    const VecT dp = p1 - p0; // Direction vector of segment S1
    const VecT dq = q1 - q0; // Direction vector of segment S2
    const VecT r = p0 - q0;

    const ScalarT a = dp.dot( dp ); // Squared length of segment S1, always nonnegative
    const ScalarT e = dq.dot( dq ); // Squared length of segment S2, always nonnegative
    const ScalarT f = dq.dot( r );

    const ScalarT c = dp.dot( r );
    const ScalarT b = dp.dot( dq );

    const ScalarT denom = a * e - b * b;

    if ( isSmall( denom ) ) // parallel
    {
        const ScalarT s0 = -c / a;
        const ScalarT s1 = ( b - c ) / a;

        s = -1; // FIXME
        t = -1; // FIXME

        if ( s0 < 0 )
        {
            s = 0;
            if ( s1 < s0 )
            {
                t = 0;
                return r.squaredNorm();
            }
            else if ( s1 < 0 )
            {
                t = 1;
                return ( p0 - q1 ).squaredNorm();
            }
        }
        else if ( s0 > 1 )
        {
            s = 1;
            if ( s1 > s0 )
            {
                t = 0;
                return ( p1 - q0 ).squaredNorm();
            }
            else if ( s1 > 1 )
            {
                t = 1;
                return ( p1 - q1 ).squaredNorm();
            }
        }

        return ( r - c * dp / a ).squaredNorm();
    }

    const ScalarT s_ = ( b * f - c * e ) / denom;
    const ScalarT t_ = ( b * s_ + f ) / e;

    s = clamp<ScalarT>( s_, 0, 1 );
    t = clamp<ScalarT>( t_, 0, 1 );

    return ( p0 + s * dp - q0 - t * dq ).squaredNorm();
}

template<typename VecT, typename ScalarT, typename InVecT>
ScalarT DistSegmentToSegment( const Eigen::MatrixBase<InVecT>& p0,
        const Eigen::MatrixBase<InVecT>& p1, const Eigen::MatrixBase<InVecT>& q0,
        const Eigen::MatrixBase<InVecT>& q1 )
{
    ScalarT s, t;
    std::sqrt( SquareDistSegmentToSegment<VecT, ScalarT, InVecT>( p0, p1, q0, q1, s, t ) );
}

template<typename VecT, typename ScalarT, typename InVecT>
ScalarT HausdorffDistSegmentToSegment( const Eigen::MatrixBase<InVecT>& p0,
        const Eigen::MatrixBase<InVecT>& p1, const Eigen::MatrixBase<InVecT>& q0,
        const Eigen::MatrixBase<InVecT>& q1 )
{
    const VecT dp = p1 - p0; // Direction vector of segment S1
    const VecT dq = q1 - q0; // Direction vector of segment S2
    const VecT r = p0 - q0;

    const ScalarT a = dp.dot( dp ); // Squared length of segment S1, always nonnegative
    const ScalarT e = dq.dot( dq ); // Squared length of segment S2, always nonnegative

    const ScalarT f = dq.dot( r );
    const ScalarT c = dp.dot( r );
    const ScalarT b = dp.dot( dq );

    const ScalarT s0 = clamp( -c / a, 0.f, 1.f );
    const ScalarT s1 = clamp( ( b - c ) / a, 0.f, 1.f );

    const ScalarT t0 = clamp( f / e, 0.f, 1.f );
    const ScalarT t1 = clamp( ( b + f ) / e, 0.f, 1.f );

    return std::sqrt(
            std::max(
                    std::max( ( p0 + s0 * dp - q0 ).squaredNorm(),
                            ( p0 + s1 * dp - q1 ).squaredNorm() ),
                    std::max( ( q0 + t0 * dq - p0 ).squaredNorm(),
                            ( q0 + t1 * dq - p1 ).squaredNorm() ) ) );

}

#endif 
