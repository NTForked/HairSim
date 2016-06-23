#include "ElasticStrandUtils.hh"

Vec3x parallelTransport( const Vec3x& u, const Vec3x& t0, const Vec3x& t1 )
{
    // Compute rotation axis (if any)
    Vec3x b = t0.cross( t1 );
    const Scalar bNorm = b.norm();
    if ( isSmall( bNorm ) ) // vectors are nearly collinear
        return u;
    b /= bNorm;

    const Vec3x& n0 = t0.cross( b ).normalized();
    const Vec3x& n1 = t1.cross( b ).normalized();

    return u.dot( t0.normalized() ) * t1.normalized() + u.dot( n0 ) * n1 + u.dot( b ) * b;
}

Vec3x normalParallelTransport( const Vec3x& u, const Vec3x& t0, const Vec3x& t1 )
{
    // This should be called only to transport an orthogonal vector
    assert( isSmall(u.dot(t0)) );

    // Compute rotation axis (if any)
    Vec3x b = t0.cross( t1 );
    const Scalar bNorm = b.norm();
    if ( isSmall( bNorm ) ) // vectors are nearly collinear
        return u;
    b /= bNorm;

    const Vec3x& n0 = t0.cross( b ).normalized();
    const Vec3x& n1 = t1.cross( b ).normalized();

    return u.dot( n0 ) * n1 + u.dot( b ) * b;
}

Vec3x orthonormalParallelTransport( const Vec3x& u, const Vec3x& t0, const Vec3x& t1 )
{
    // This should be called only to transport an orthogonal vector
    assert( isSmall(u.dot(t0)) );

    Vec3x b = t0.cross( t1 );
    const Scalar bNorm = b.norm();
    if ( isSmall( bNorm ) ) // vectors are nearly collinear
        return u;
    b /= bNorm;

    const Vec3x& n0 = t0.cross( b );
    const Vec3x& n1 = t1.cross( b );

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
