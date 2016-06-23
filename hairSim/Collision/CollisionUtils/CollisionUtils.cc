#include "CollisionUtils.hh"
#include "OrientedBoundingBox.hh"

#include "../Core/ElasticStrand.hh"
#include "../Utils/Distances.hh"

#include "../Dynamic/StrandDynamicTraits.hh"
#include "../Dynamic/Config.hh"

#define COS_PARALLEL_ENOUGH 0.86602540378      // cos( Pi/6 )

// Adapted from some code on Robert Bridson's website, I believe

void addUnique( std::vector<double>& a, double e )
{
    for ( unsigned int i = 0; i < a.size(); ++i )
        if ( a[i] == e )
            return;
    a.push_back( e );
}

void addUnique( double *a, unsigned &a_size, double e )
{
    for ( unsigned int i = 0; i < a_size; ++i )
        if ( a[i] == e )
            return;
    a[a_size++] = e;
}

inline void compare_and_swap( double& a, double&b )
{
    if ( a > b )
        std::swap( a, b );
}

//! Sorting networks for a_size <= 4
void sort( double *a, unsigned a_size )
{
    switch ( a_size )
    {
    case 4:
        compare_and_swap( a[0], a[2] );
        compare_and_swap( a[1], a[3] );
        compare_and_swap( a[0], a[1] );
        compare_and_swap( a[2], a[3] );
        compare_and_swap( a[1], a[2] );
        break;
    case 3:
        compare_and_swap( a[0], a[1] );
        compare_and_swap( a[0], a[2] );
        compare_and_swap( a[1], a[2] );
        break;
    case 2:
        compare_and_swap( a[0], a[1] );
        break;
    default:
        break;
    }
}

double triple( const Vec3x& a, const Vec3x& b, const Vec3x& c )
{
    return a[0] * ( b[1] * c[2] - b[2] * c[1] ) + a[1] * ( b[2] * c[0] - b[0] * c[2] )
            + a[2] * ( b[0] * c[1] - b[1] * c[0] );
}

double signed_volume( const Vec3x& x0, const Vec3x& x1, const Vec3x& x2, const Vec3x& x3 )
{
    // Equivalent to triple(x1-x0, x2-x0, x3-x0), six times the signed volume of the tetrahedron.
    // But, for robustness, we want the result (up to sign) to be independent of the ordering.
    // And want it as accurate as possible...
    // But all that stuff is hard, so let's just use the common assumption that all coordinates are >0,
    // and do something reasonably accurate in fp.

    // This formula does almost four times too much multiplication, but if the coordinates are non-negative
    // it suffers in a minimal way from cancellation error.
    return ( x0[0] * ( x1[1] * x3[2] + x3[1] * x2[2] + x2[1] * x1[2] )
            + x1[0] * ( x2[1] * x3[2] + x3[1] * x0[2] + x0[1] * x2[2] )
            + x2[0] * ( x3[1] * x1[2] + x1[1] * x0[2] + x0[1] * x3[2] )
            + x3[0] * ( x1[1] * x2[2] + x2[1] * x0[2] + x0[1] * x1[2] ) )

            - ( x0[0] * ( x2[1] * x3[2] + x3[1] * x1[2] + x1[1] * x2[2] )
                    + x1[0] * ( x3[1] * x2[2] + x2[1] * x0[2] + x0[1] * x3[2] )
                    + x2[0] * ( x1[1] * x3[2] + x3[1] * x0[2] + x0[1] * x1[2] )
                    + x3[0] * ( x2[1] * x1[2] + x1[1] * x0[2] + x0[1] * x2[2] ) );
}

// All roots returned in interval [0,1]. Assumed geometry followed a linear
// trajectory between x and xnew. 
void getCoplanarityTimes( const Vec3x& x0, const Vec3x& x1, const Vec3x& x2, const Vec3x& x3,
        const Vec3x& xnew0, const Vec3x& xnew1, const Vec3x& xnew2, const Vec3x& xnew3,
        double* times, double* errors, unsigned &num_times )
{
    const double tol = 1e-8;
    num_times = 0;

    if( thickness ){ 
        times[num_times++] = 0.; // TODO: unclear if these should be here still if thickness is off.
        times[num_times++] = 1.; // are there still boundary conditions to check?
    }

    // cubic coefficients, A*t^3+B*t^2+C*t+D (for t in [0,1])
    const Vec3x x03 = x0 - x3;
    const Vec3x x13 = x1 - x3;
    const Vec3x x23 = x2 - x3;
    const Vec3x v03 = ( xnew0 - xnew3 ) - x03;
    const Vec3x v13 = ( xnew1 - xnew3 ) - x13;
    const Vec3x v23 = ( xnew2 - xnew3 ) - x23;

    double A = triple( v03, v13, v23 );
    double B = triple( x03, v13, v23 ) + triple( v03, x13, v23 ) + triple( v03, v13, x23 );
    double C = triple( x03, x13, v23 ) + triple( x03, v13, x23 ) + triple( v03, x13, x23 );
    double D = triple( x03, x13, x23 );

    const double convergence_tol = tol
            * ( std::fabs( A ) + std::fabs( B ) + std::fabs( C ) + std::fabs( D ) );

    // find intervals to check, or just solve it if it reduces to a quadratic =============================
    double interval_times[4];
    unsigned interval_times_size = 0;

    double discriminant = B * B - 3 * A * C; // of derivative of cubic, 3*A*t^2+2*B*t+C, divided by 4 for convenience
    if ( discriminant <= 0 )
    { // monotone cubic: only one root in [0,1] possible
      // so we just
        interval_times[0] = 0;
        interval_times[1] = 1;
        interval_times_size = 2;
    }
    else
    { // positive discriminant, B!=0
        if ( A == 0 )
        { // the cubic is just a quadratic, B*t^2+C*t+D ========================================
            discriminant = C * C - 4 * B * D; // of the quadratic
            if ( discriminant <= 0 )
            {
                double t = -C / ( 2 * B );
                if ( t >= -tol && t <= 1 + tol )
                {
                    t = clamp( t, 0., 1. );
                    double val = std::fabs(
                            signed_volume( ( 1 - t ) * x0 + t * xnew0, ( 1 - t ) * x1 + t * xnew1,
                                    ( 1 - t ) * x2 + t * xnew2, ( 1 - t ) * x3 + t * xnew3 ) );
                    if ( val < convergence_tol )
                    {
                        times[num_times++] = t;
                    }
                }
            }
            else
            { // two separate real roots
                double t0, t1;
                if ( C > 0 )
                    t0 = ( -C - std::sqrt( discriminant ) ) / ( 2 * B );
                else
                    t0 = ( -C + std::sqrt( discriminant ) ) / ( 2 * B );
                t1 = D / ( B * t0 );
                if ( t1 < t0 )
                    std::swap( t0, t1 );
                if ( t0 >= -tol && t0 <= 1 + tol )
                {
                    times[num_times++] = clamp( t0, 0., 1. );
                }
                if ( t1 >= -tol && t1 <= 1 + tol )
                {
                    addUnique( times, num_times, clamp( t1, 0., 1. ) );
                }
            }

            if ( errors )
            {
                for ( unsigned i = 0; i < num_times; ++i )
                {
                    double ti = times[i];
                    double val = std::fabs(
                            signed_volume( ( 1 - ti ) * x0 + ti * xnew0,
                                    ( 1 - ti ) * x1 + ti * xnew1, ( 1 - ti ) * x2 + ti * xnew2,
                                    ( 1 - ti ) * x3 + ti * xnew3 ) );
                    errors[i] = val;
                }
            }

            return;
        }
        else
        { // cubic is not monotone: divide up [0,1] accordingly =====================================
            double t0, t1;
            if ( B > 0 )
                t0 = ( -B - std::sqrt( discriminant ) ) / ( 3 * A );
            else
                t0 = ( -B + std::sqrt( discriminant ) ) / ( 3 * A );
            t1 = C / ( 3 * A * t0 );
            if ( t1 < t0 )
                std::swap( t0, t1 );

            interval_times[interval_times_size++] = 0;
            if ( t0 > 0 && t0 < 1 )
                interval_times[interval_times_size++] = t0;
            if ( t1 > 0 && t1 < 1 )
                interval_times[interval_times_size++] = t1;

            interval_times[interval_times_size++] = 1;
        }
    }

    // look for roots in indicated intervals ==============================================================
    // evaluate coplanarity more accurately at each endpoint of the intervals
    double interval_values[interval_times_size];
    for ( unsigned int i = 0; i < interval_times_size; ++i )
    {
        double t = interval_times[i];
        interval_values[i] = signed_volume( ( 1 - t ) * x0 + t * xnew0, ( 1 - t ) * x1 + t * xnew1,
                ( 1 - t ) * x2 + t * xnew2, ( 1 - t ) * x3 + t * xnew3 );
    }
    // first look for interval endpoints that are close enough to zero, without a sign change
    for ( unsigned int i = 0; i < interval_times_size; ++i )
    {
        if ( interval_values[i] == 0 )
        {
            times[num_times++] = interval_times[i];
        }
        else if ( std::fabs( interval_values[i] ) < convergence_tol )
        {
            if ( ( i == 0 || ( interval_values[i - 1] >= 0 && interval_values[i] >= 0 )
                    || ( interval_values[i - 1] <= 0 && interval_values[i] <= 0 ) )
                    && ( i == interval_times_size - 1
                            || ( interval_values[i + 1] >= 0 && interval_values[i] >= 0 )
                            || ( interval_values[i + 1] <= 0 && interval_values[i] <= 0 ) ) )
            {
                times[num_times++] = interval_times[i];
            }
        }
    }
    // and then search in intervals with a sign change
    for ( unsigned int i = 1; i < interval_times_size; ++i )
    {
        double tlo = interval_times[i - 1], thi = interval_times[i], tmid;
        double vlo = interval_values[i - 1], vhi = interval_values[i], vmid;
        if ( ( vlo < 0 && vhi > 0 ) || ( vlo > 0 && vhi < 0 ) )
        {
            // start off with secant approximation (in case the cubic is actually linear)
            double alpha = vhi / ( vhi - vlo );
            tmid = alpha * tlo + ( 1 - alpha ) * thi;
            for ( int iteration = 0; iteration < 50; ++iteration )
            {
                vmid = signed_volume( ( 1 - tmid ) * x0 + tmid * xnew0,
                        ( 1 - tmid ) * x1 + tmid * xnew1, ( 1 - tmid ) * x2 + tmid * xnew2,
                        ( 1 - tmid ) * x3 + tmid * xnew3 );
                if ( std::fabs( vmid ) < 1e-2 * convergence_tol )
                    break;
                if ( ( vlo < 0 && vmid > 0 ) || ( vlo > 0 && vmid < 0 ) )
                { // if sign change between lo and mid
                    thi = tmid;
                    vhi = vmid;
                }
                else
                { // otherwise sign change between hi and mid
                    tlo = tmid;
                    vlo = vmid;
                }
                if ( iteration % 2 )
                    alpha = 0.5; // sometimes go with bisection to guarantee we make progress
                else
                    alpha = vhi / ( vhi - vlo ); // other times go with secant to hopefully get there fast
                tmid = alpha * tlo + ( 1 - alpha ) * thi;
            }
            times[num_times++] = tmid;
        }
    }
    sort( times, num_times );

    if ( errors )
    {
        for ( unsigned i = 0; i < num_times; ++i )
        {
            double ti = times[i];
            double val = std::fabs(
                    signed_volume( ( 1 - ti ) * x0 + ti * xnew0, ( 1 - ti ) * x1 + ti * xnew1,
                            ( 1 - ti ) * x2 + ti * xnew2, ( 1 - ti ) * x3 + ti * xnew3 ) );
            errors[i] = val;
        }
    }
}

void getIntersectionPoint( const Vec3x& x_edge_0, const Vec3x& x_edge_1, const Vec3x& x_face_0,
        const Vec3x& x_face_1, const Vec3x& x_face_2, double* times, double* errors,
        unsigned &num_times )
{
    const double tol = 1e-12;
    num_times = 0;

    Vec3x x03 = x_edge_0 - x_face_2;
    Vec3x x13 = x_face_0 - x_face_2;
    Vec3x x23 = x_face_1 - x_face_2;
    Vec3x v03 = x_edge_1 - x_face_2 - x03;

    double C = triple( v03, x13, x23 );
    double D = triple( x03, x13, x23 );

    const double convergence_tol = tol
            * ( std::fabs( 0 ) + std::fabs( 0 ) + std::fabs( C ) + std::fabs( D ) );

    // find intervals to check, or just solve it if it reduces to a quadratic =============================
    const unsigned interval_times_size = 2;
    double interval_times[interval_times_size] = { 0., 1. };

    // look for roots in indicated intervals ==============================================================
    // evaluate coplanarity more accurately at each endpoint of the intervals
    double interval_values[interval_times_size];

    for ( unsigned int i = 0; i < interval_times_size; ++i )
    {
        double t = interval_times[i];
        interval_values[i] = signed_volume( ( 1 - t ) * x_edge_0 + t * x_edge_1, x_face_0, x_face_1,
                x_face_2 );
    }
    // first look for interval endpoints that are close enough to zero, without a sign change
    for ( unsigned int i = 0; i < interval_times_size; ++i )
    {
        if ( interval_values[i] == 0 )
        {
            times[num_times++] = interval_times[i];
        }
        else if ( std::fabs( interval_values[i] ) < convergence_tol )
        {
            if ( ( i == 0 || ( interval_values[i - 1] >= 0 && interval_values[i] >= 0 )
                    || ( interval_values[i - 1] <= 0 && interval_values[i] <= 0 ) )
                    && ( i == interval_times_size - 1
                            || ( interval_values[i + 1] >= 0 && interval_values[i] >= 0 )
                            || ( interval_values[i + 1] <= 0 && interval_values[i] <= 0 ) ) )
            {
                times[num_times++] = interval_times[i];
            }
        }
    }
    // and then search in intervals with a sign change
    for ( unsigned int i = 1; i < interval_times_size; ++i )
    {
        double tlo = interval_times[i - 1], thi = interval_times[i], tmid;
        double vlo = interval_values[i - 1], vhi = interval_values[i], vmid;
        if ( ( vlo < 0 && vhi > 0 ) || ( vlo > 0 && vhi < 0 ) )
        {
            // start off with secant approximation (in case the cubic is actually linear)
            double alpha = vhi / ( vhi - vlo );
            tmid = alpha * tlo + ( 1 - alpha ) * thi;
            for ( int iteration = 0; iteration < 50; ++iteration )
            {
                vmid = signed_volume( ( 1 - tmid ) * x_edge_0 + tmid * x_edge_1,
                        ( 1 - tmid ) * x_face_0 + tmid * x_face_0,
                        ( 1 - tmid ) * x_face_1 + tmid * x_face_1,
                        ( 1 - tmid ) * x_face_2 + tmid * x_face_2 );
                if ( std::fabs( vmid ) < 1e-2 * convergence_tol )
                    break;
                if ( ( vlo < 0 && vmid > 0 ) || ( vlo > 0 && vmid < 0 ) )
                { // if sign change between lo and mid
                    thi = tmid;
                    vhi = vmid;
                }
                else
                { // otherwise sign change between hi and mid
                    tlo = tmid;
                    vlo = vmid;
                }
                if ( iteration % 2 )
                    alpha = 0.5; // sometimes go with bisection to guarantee we make progress
                else
                    alpha = vhi / ( vhi - vlo ); // other times go with secant to hopefully get there fast
                tmid = alpha * tlo + ( 1 - alpha ) * thi;
            }
            times[num_times++] = tmid;
        }
    }
    sort( times, num_times );

    if ( errors )
    {
        for ( unsigned i = 0; i < num_times; ++i )
        {
            double ti = times[i];
            double val = std::fabs(
                    signed_volume( ( 1 - ti ) * x_edge_0 + ti * x_edge_1, x_face_0, x_face_1,
                            x_face_2 ) );
            errors[i] = val;
        }
    }
}

void buildFrame( const Vec3x& n_hat, Vec3x& t1, Vec3x& t2 )
{
    // build collision frame basis (n_hat, t1, t2)
    Vec3x tmp( n_hat( 0 ), -n_hat( 2 ), n_hat( 1 ) );
    if ( n_hat( 2 ) * n_hat( 2 ) < 1e-12 && n_hat( 1 ) * n_hat( 1 ) < 1e-12 )
        tmp << -n_hat( 1 ), n_hat( 0 ), n_hat( 2 );
    t1 = tmp.cross( n_hat );

    if ( t1.norm() < 1e-8 )
    {
        ErrorStream( g_log, "" ) << "tangent degenerate: tmp =" << tmp << ", n = " << n_hat
                << ", t1 =" << t1;
        return;
    }

    t1.normalize();
    t2 = t1.cross( n_hat );
//    t2 = t2 / t2.norm();
    assert( isApproxUnit( t2 ) );
}

bool analyseRoughRodRodCollision( const ElasticStrand* sP, const ElasticStrand* sQ, const int iP,
        const int iQ, Vec3x &depl, Scalar &s, Scalar &t, Scalar &d )
{
    const CollisionParameters &cpP = sP->collisionParameters();
    const CollisionParameters &cpQ = sQ->collisionParameters();

    const Vec3x &P0 = sP->getVertex( iP );
    const Vec3x &P1 = sP->getVertex( iP + 1 );
    const Vec3x &Q0 = sQ->getVertex( iQ );
    const Vec3x &Q1 = sQ->getVertex( iQ + 1 );

    const Scalar BCRad = cpP.collisionRadius( iP ) + cpQ.collisionRadius( iQ );

    Scalar sqDist = SquareDistSegmentToSegment<Vec3x, Scalar, Vec3x>( P0, P1, Q0, Q1, s, t );
    // Required to determnisticity -- are x87 registers sometimes used ?
    s = ( float ) s;
    t = ( float ) t;

    if ( sqDist > BCRad * BCRad || s < 0 || t < 0 )
        return false; // see FIXME in DistSegmentToSegment

    Vec3x PC = ( ( 1. - s ) * P0 + s * P1 );
    Vec3x QC = ( ( 1. - t ) * Q0 + t * Q1 );
    depl = ( PC - QC );

    const Scalar n2depl = depl.squaredNorm();
    if ( isSmall( n2depl ) )
        return false;

    depl /= std::sqrt( n2depl );

    if ( depl.dot( ( P1 - P0 ).normalized() ) > COS_PARALLEL_ENOUGH
            || depl.dot( ( Q1 - Q0 ).normalized() ) < -COS_PARALLEL_ENOUGH )
        return false;

    d = std::sqrt( sqDist );

    return true;
}

static const double SQ_TOLERANCE = 1e-12;

void normalFunc(
    const Vec3x &p1_old, const Vec3x &q1_old, 
    const Vec3x &p2_old, const Vec3x &q2_old, 
    const Vec3x &p1_new, const Vec3x &q1_new, 
    const Vec3x &p2_new, const Vec3x &q2_new,
    const Scalar alpha1, const Scalar alpha2, Vec3x& normal ) 
{
    Vec3x e1_new = q1_new - p1_new;
    Vec3x e2_new = q2_new - p2_new;
    assert(SQ_TOLERANCE <= e1_new.squaredNorm());
    assert(SQ_TOLERANCE <= e2_new.squaredNorm());
    e1_new.normalize();
    e2_new.normalize();

    // if e1 and e2 are not parallel, use their cross product
    Vec3x n  = e1_new.cross(e2_new);
    if (n.squaredNorm() < SQ_TOLERANCE) {
        // e1 and e2 are parallel, so find a vector perp. to e1
        // try starting with the vector connecting the intersection points
        // as observed at start of time step
        n = ((1 - alpha2) * p2_old + alpha2 * q2_old) - 
            ((1 - alpha1) * p1_old + alpha1 * q1_old);

        // project out component along e1_new
        n -= (n.dot(e1_new)) * e1_new;

        // if that fails, try other vectors
        if (n.squaredNorm() < SQ_TOLERANCE) {
            n = Vec3x(1,0,0);
            n -= (n.dot(e1_new)) * e1_new;

            if (n.squaredNorm() < SQ_TOLERANCE ){
                n = Vec3x(0,1,0);
                n -= (n.dot(e1_new)) * e1_new;
            }
        }
    }            
    assert(SQ_TOLERANCE <= n.squaredNorm());
    n.normalize();    

    // let x0 and x1 be the full (12d) configuration at the start of time step
    // and at the instant of collision, respectively.
    //
    // Let's adopt the convention that the 12d space is ordered (p1,q1,p2,q2)
    //
    // Compute the normal velocity nvel = (x1-x0).dot(n),
    // where n12 is the 12d (configurational) normal at the time of collision,
    // i.e., n12 = (-(1-alpha1)*n, -alpha1*n, (1-alpha2)*n, alpha2*n)

    typedef Eigen::Matrix<Scalar, 12, 1>  Vec12;
    
    Vec12 x_new;
    x_new <<
        p1_new(0), p1_new(1), p1_new(2),
        q1_new(0), q1_new(1), q1_new(2),
        p2_new(0), p2_new(1), p2_new(2),
        q2_new(0), q2_new(1), q2_new(2);

    Vec12 x_old;
    x_old <<
        p1_old(0), p1_old(1), p1_old(2),
        q1_old(0), q1_old(1), q1_old(2),
        p2_old(0), p2_old(1), p2_old(2),
        q2_old(0), q2_old(1), q2_old(2);

    Vec12 n12;
    n12 <<
        -(1-alpha1)*n(0), -(1-alpha1)*n(1), -(1-alpha1)*n(2),
        -(  alpha1)*n(0), -(  alpha1)*n(1), -(  alpha1)*n(2),
        +(1-alpha2)*n(0), +(1-alpha2)*n(1), +(1-alpha2)*n(2),
        +(  alpha2)*n(0), +(  alpha2)*n(1), +(  alpha2)*n(2);
    n12.normalize();

    const Vec12 x_col = x_new - x_new.dot(n12)*n12;
    
    const Scalar normal_distance = (x_old - x_col).dot(n12);

    if (normal_distance < 0.) {
        n = -n;
    }

    normal = n;
}

bool compareCT( const CollisionBase* ct1, const CollisionBase* ct2 )
{
    {
        const VertexFaceCollision* vf1 = dynamic_cast<const VertexFaceCollision*>( ct1 );
        const VertexFaceCollision* vf2 = dynamic_cast<const VertexFaceCollision*>( ct2 );

        if( vf1 != NULL && vf2 == NULL ){
            return true;
        }
        if( vf1 == NULL && vf2 != NULL ){
            return false;
        }
        if( vf1 != NULL && vf2 != NULL ){
            return compare( vf1, vf2 );
        }
    }

    {
        const EdgeFaceCollision* ef1 = dynamic_cast<const EdgeFaceCollision*>( ct1 );
        const EdgeFaceCollision* ef2 = dynamic_cast<const EdgeFaceCollision*>( ct2 );

        if( ef1 != NULL && ef2 == NULL ){
            return true;
        }
        if( ef1 == NULL && ef2 != NULL ){
            return false;
        }
        if( ef1 != NULL && ef2 != NULL ){
            return compare( ef1, ef2 );
        }
    }

    {
        const EdgeEdgeCollision* ee1 = dynamic_cast<const EdgeEdgeCollision*>( ct1 );
        const EdgeEdgeCollision* ee2 = dynamic_cast<const EdgeEdgeCollision*>( ct2 );

        if( ee1 != NULL && ee2 == NULL ){
            return true;
        }
        if( ee1 == NULL && ee2 != NULL ){
            return false;
        }
        if( ee1 != NULL && ee2 != NULL ){
            return compare( ee1, ee2 );
        }
    }
    std::cerr << "We should never arrive here!\n";
    return false;
}