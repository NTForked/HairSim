#include "ElasticStrandUtils.hh"

Vec3 parallelTransport( const Vec3& u, const Vec3& t0, const Vec3& t1 )
{
    // Compute rotation axis (if any)
    Vec3 b = t0.cross( t1 );
    const Scalar bNorm = b.norm();
    if ( isSmall( bNorm ) ) // vectors are nearly collinear
        return u;
    b /= bNorm;

    const Vec3& n0 = t0.cross( b ).normalized();
    const Vec3& n1 = t1.cross( b ).normalized();

    return u.dot( t0.normalized() ) * t1.normalized() + u.dot( n0 ) * n1 + u.dot( b ) * b;
}

Vec3 normalParallelTransport( const Vec3& u, const Vec3& t0, const Vec3& t1 )
{
    // This should be called only to transport an orthogonal vector
    assert( isSmall(u.dot(t0)) );

    // Compute rotation axis (if any)
    Vec3 b = t0.cross( t1 );
    const Scalar bNorm = b.norm();
    if ( isSmall( bNorm ) ) // vectors are nearly collinear
        return u;
    b /= bNorm;

    const Vec3& n0 = t0.cross( b ).normalized();
    const Vec3& n1 = t1.cross( b ).normalized();

    return u.dot( n0 ) * n1 + u.dot( b ) * b;
}

Vec3 orthonormalParallelTransport( const Vec3& u, const Vec3& t0, const Vec3& t1 )
{
    // This should be called only to transport an orthogonal vector
    assert( isSmall(u.dot(t0)) );

    Vec3 b = t0.cross( t1 );
    const Scalar bNorm = b.norm();
    if ( isSmall( bNorm ) ) // vectors are nearly collinear
        return u;
    b /= bNorm;

    const Vec3& n0 = t0.cross( b );
    const Vec3& n1 = t1.cross( b );

    return u.dot( n0 ) * n1 + u.dot( b ) * b;
}

bool containsNans( const VecXx &dofs )
{
    bool panic = false;
    for ( int i = 0; i < dofs.size(); ++i )
    {
        if ( std::isnan( dofs[i] ) )
        {
            panic = true;
            break;
        }
    }

    return panic;
}
