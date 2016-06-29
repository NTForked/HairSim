#include "StrandStrandForce.hh"

#include "../Static/CollisionSet.hh"
#include "../Static/StrandStaticTraits.hh"
#include "../Core/ElasticStrand.hh"
#include "../Core/BandMatrix.hh"

//#define USE_DISTANCE_FOR_PENALTY_FORCES

StrandStrandForce::StrandStrandForce( const ElasticStrand* strand1, int edge1,
        const ElasticStrand* strand2, int edge2, double thickness, double stiffness,
        double hardness, double friction ) :
        RuntimeCollisionForceBase(), 
        m_thickness( thickness ), 
        m_stiffness( stiffness ), 
        m_hardness( hardness ), 
        m_friction( friction ), 
        m_colliding( false ), 
        m_strandP_col_initialized( false ), 
        m_strandQ_col_initialized( false )
{
    if ( strand1 == strand2 )
    {
        if ( edge1 < edge2 )
        {
            m_edgeP = edge1;
            m_edgeQ = edge2;
        }
        else
        {
            m_edgeP = edge2;
            m_edgeQ = edge1;
        }
    }
    else if ( strand1 > strand2 )
    {
        m_strandP = strand1;
        m_strandQ = strand2;
        m_edgeP = edge1;
        m_edgeQ = edge2;
    }
    else
    {
        m_strandP = strand2;
        m_strandQ = strand1;
        m_edgeP = edge2;
        m_edgeQ = edge1;
    }

    // Compute the initial sign of the tetrahedron's volume

    findTheVertices( strand1->getCurrentState() );
    const Scalar det = tetraVolume();

    if ( det > 0 )
        m_initialSign = +1.0;
    else if ( det < 0 )
        m_initialSign = -1.0;
    else
        m_initialSign = 0.0;

}

StrandStrandForce::StrandStrandForce( const ElasticStrand* strand1, int edge1,
        const ElasticStrand* strand2, int edge2, double thickness, double stiffness,
        double hardness, double friction, double initialSign ) :
        RuntimeCollisionForceBase(),
        m_initialSign( initialSign ),  
        m_thickness( thickness ), 
        m_stiffness( stiffness ), 
        m_hardness( hardness ), 
        m_friction( friction ), 
        m_colliding( false ), 
        m_strandP_col_initialized( false ), 
        m_strandQ_col_initialized( false )
{
    if ( strand1 == strand2 )
    {
        if ( edge1 < edge2 )
        {
            m_edgeP = edge1;
            m_edgeQ = edge2;
        }
        else
        {
            m_edgeP = edge2;
            m_edgeQ = edge1;
        }
    }
    else if ( strand1 > strand2 )
    {
        m_strandP = strand1;
        m_strandQ = strand2;
        m_edgeP = edge1;
        m_edgeQ = edge2;
    }
    else
    {
        m_strandP = strand2;
        m_strandQ = strand1;
        m_edgeP = edge2;
        m_edgeQ = edge1;
    }
}

StrandStrandForce::~StrandStrandForce()
{
    // TODO Auto-generated destructor stub
    // std::cerr << "Strand/strand destructed" << std::endl ;
}

std::string StrandStrandForce::getName() const
{
    return "Strand/strand";
}

Scalar StrandStrandForce::tetraVolume() const
{
//    const Vec3x& p0 = m_strandP->getVertex( m_edgeP );
//    const Vec3x& p1 = m_strandP->getVertex( m_edgeP + 1 );
//    const Vec3x& q0 = m_strandQ->getVertex( m_edgeQ );
//    const Vec3x& q1 = m_strandQ->getVertex( m_edgeQ + 1 );

    Mat3x tetra;
    tetra.block < 3, 1 > ( 0, 0 ) = m_p1 - m_p0;
    tetra.block < 3, 1 > ( 0, 1 ) = m_q0 - m_p0;
    tetra.block < 3, 1 > ( 0, 2 ) = m_q1 - m_p0;

//    std::cout << "Volume: " << tetra.determinant() << '\n';

    return tetra.determinant(); // That's actually 6 x the tetrahedron volume.
}

Vec3x StrandStrandForce::gradP0TetraVolume() const
{
//    const Vec3x& p1 = m_strandP->getVertex( m_edgeP + 1 );
//    const Vec3x& q0 = m_strandQ->getVertex( m_edgeQ );
//    const Vec3x& q1 = m_strandQ->getVertex( m_edgeQ + 1 );
    const Vec3x& v0 = m_p1 - m_q0;
    const Vec3x& v1 = m_q1 - m_p1;

    Vec3x grad;
    grad[0] = v0[1] * v1[2] - v0[2] * v1[1];
    grad[1] = v0[2] * v1[0] - v0[0] * v1[2];
    grad[2] = v0[0] * v1[1] - v0[1] * v1[0];

    return grad;
}

Vec3x StrandStrandForce::gradP1TetraVolume() const
{
//    const Vec3x& p0 = m_strandP->getVertex( m_edgeP );
//    const Vec3x& q0 = m_strandQ->getVertex( m_edgeQ );
//    const Vec3x& q1 = m_strandQ->getVertex( m_edgeQ + 1 );
    const Vec3x& v0 = m_p0 - m_q0;
    const Vec3x& v1 = m_q1 - m_p0;

    Vec3x grad;
    grad[0] = v0[2] * v1[1] - v0[1] * v1[2];
    grad[1] = v0[0] * v1[2] - v0[2] * v1[0];
    grad[2] = v0[1] * v1[0] - v0[0] * v1[1];

    return grad;
}

Vec3x StrandStrandForce::gradQ1TetraVolume() const
{
//    const Vec3x& p0 = m_strandP->getVertex( m_edgeP );
//    const Vec3x& p1 = m_strandP->getVertex( m_edgeP + 1 );
//    const Vec3x& q0 = m_strandQ->getVertex( m_edgeQ );
    const Vec3x& v0 = m_q0 - m_p0;
    const Vec3x& v1 = m_p1 - m_q0;

    Vec3x grad;
    grad[0] = v0[2] * v1[1] - v0[1] * v1[2];
    grad[1] = v0[0] * v1[2] - v0[2] * v1[0];
    grad[2] = v0[1] * v1[0] - v0[0] * v1[1];

    return grad;
}

Mat3x StrandStrandForce::hessP0P1TetraVolume() const
{
//    const Vec3x& q0 = m_strandQ->getVertex( m_edgeQ );
//    const Vec3x& q1 = m_strandQ->getVertex( m_edgeQ + 1 );

    Mat3x hess;
    hess( 0, 0 ) = hess( 1, 1 ) = hess( 2, 2 ) = 0;
    hess( 0, 1 ) = m_q0[2] - m_q1[2];
    hess( 1, 0 ) = -hess( 0, 1 );
    hess( 0, 2 ) = m_q1[1] - m_q0[1];
    hess( 2, 0 ) = -hess( 0, 2 );
    hess( 1, 2 ) = m_q0[0] - m_q1[0];
    hess( 2, 1 ) = -hess( 1, 2 );

    return hess;
}

Mat3x StrandStrandForce::hessQ0Q1TetraVolume() const
{
//    const Vec3x& p0 = m_strandQ->getVertex( m_edgeP );
//    const Vec3x& p1 = m_strandQ->getVertex( m_edgeP + 1 );

    Mat3x hess;
    hess( 0, 0 ) = hess( 1, 1 ) = hess( 2, 2 ) = 0;
    hess( 0, 1 ) = m_p0[2] - m_p1[2];
    hess( 1, 0 ) = -hess( 0, 1 );
    hess( 0, 2 ) = m_p1[1] - m_p0[1];
    hess( 2, 0 ) = -hess( 0, 2 );
    hess( 1, 2 ) = m_p0[0] - m_p1[0];
    hess( 2, 1 ) = -hess( 1, 2 );

    return hess;
}

void StrandStrandForce::findTheVertices( const StrandState& geometry ) const
{
    const StrandState* pGeo = &geometry;
    if ( pGeo == &( m_strandP->getFutureState() ) )
    {
        m_p0 = m_strandP->getFutureState().getVertex( m_edgeP );
        m_p1 = m_strandP->getFutureState().getVertex( m_edgeP + 1 );
        m_q0 = m_strandQ->getVertex( m_edgeQ );
        m_q1 = m_strandQ->getVertex( m_edgeQ + 1 );
    }
    else if ( &geometry == &( m_strandQ->getFutureState() ) )
    {
        m_p0 = m_strandP->getVertex( m_edgeP );
        m_p1 = m_strandP->getVertex( m_edgeP + 1 );
        m_q0 = m_strandQ->getFutureState().getVertex( m_edgeQ );
        m_q1 = m_strandQ->getFutureState().getVertex( m_edgeQ + 1 );
    }
    else
    {
        m_p0 = m_strandP->getVertex( m_edgeP );
        m_p1 = m_strandP->getVertex( m_edgeP + 1 );
        m_q0 = m_strandQ->getVertex( m_edgeQ );
        m_q1 = m_strandQ->getVertex( m_edgeQ + 1 );
    }
}

Vec3x StrandStrandForce::gradQ0TetraVolume() const
{
//    const Vec3x& p0 = m_strandP->getVertex( m_edgeP );
//    const Vec3x& p1 = m_strandP->getVertex( m_edgeP + 1 );
//    const Vec3x& q1 = m_strandQ->getVertex( m_edgeQ + 1 );
    const Vec3x& v0 = m_q1 - m_p0;
    const Vec3x& v1 = m_p1 - m_q1;

    Vec3x grad;
    grad[0] = v0[1] * v1[2] - v0[2] * v1[1];
    grad[1] = v0[2] * v1[0] - v0[0] * v1[2];
    grad[2] = v0[0] * v1[1] - v0[1] * v1[0];

    return grad;
}

void StrandStrandForce::accumulateEFJ( StrandState& geometry, const ElasticStrand& strand )
{
    accumulateEF( geometry, strand );
    accumulateJ( geometry, strand );
}

void StrandStrandForce::accumulateEF( StrandState& geometry, const ElasticStrand& strand )
{
    ClientNumber clientStrand;
    if ( &strand == m_strandP )
        clientStrand = P;
    else if ( &strand == m_strandQ )
        clientStrand = Q;
    else
        return; // Nothing to do if we are trying to apply this force to another strand. This should never happen...

    findTheVertices( geometry );

    accumulateEF( clientStrand, geometry.m_totalForce, geometry.m_totalEnergy );

    registerCollision( strand, geometry );

//    std::cout << "Client: " << clientStrand << " energy: " << localE << " force: " << localF
//            << '\n';
}

void StrandStrandForce::accumulateJ( StrandState& geometry, const ElasticStrand& strand )
{

    ClientNumber clientStrand;
    if ( &strand == m_strandP )
        clientStrand = P;
    else if ( &strand == m_strandQ )
        clientStrand = Q;
    else
        return; // Nothing to do if we are trying to apply this force to another strand. This should never happen...

    findTheVertices( geometry );

#ifdef USE_DISTANCE_FOR_PENALTY_FORCES
    Scalar s, t;
    GetClosestPoints( s, t );

    if ( s == 0. || s == 1. || t == 0. || t == 1. )
        return;

    const Vec3x &ps = ( 1. - s ) * m_p0 + s * m_p1;
    const Vec3x &qt = ( 1. - t ) * m_q0 + t * m_q1;

    Vec3x edge, otherEdge;
    if ( clientStrand == P )
    {
        edge = ( m_p1 - m_p0 );
        otherEdge = ( m_q1 - m_q0 ).normalized();
    }
    else
    {
        edge = ( m_q1 - m_q0 );
        otherEdge = ( m_p1 - m_p0 ).normalized();
    }
    const Scalar edgeLen = edge.norm();

    Mat3x cross = -m_initialSign * crossMat( otherEdge );
    if ( clientStrand == Q )
    {
        cross *= -1;
    }

    const Vec3x grad = ( cross * edge ).normalized();
    const double rad = ( ps - qt ).dot( grad );
    const double depl = m_thickness - rad;

    if ( depl > 0. )
    {
        // X = P|Q
        const Mat3x dGrad_dX = ( Mat3x::Identity() - grad * grad.transpose() ) * cross;
        const Mat3x dGrad_dX0 = -dGrad_dX;
        const Mat3x dGrad_dX1 = dGrad_dX;

        Mat3x dGradGrad_dGrad[3];
        for ( unsigned i = 0; i < 3; ++i )
        {
            dGradGrad_dGrad[i].setZero();
            dGradGrad_dGrad[i].col( i ) += grad;
            dGradGrad_dGrad[i].row( i ) += grad;
        }

        Mat3x dGrad_dX_dXi[3];
        for ( unsigned i = 0; i < 3; ++i )
        {
            Mat3x dGradGrad_dXi;
            dGradGrad_dXi.setZero();

            for ( unsigned j = 0; j < 3; ++j )
            {
                dGradGrad_dXi += dGradGrad_dGrad[j] * dGrad_dX( j, i );
            }

            dGrad_dX_dXi[i] = dGradGrad_dXi * cross;
        }

        // alpha = d(Ps - Qt) / dX
        Scalar dPQ_dX0 = clientStrand == P ? ( 1 - s) : ( t - 1);
        Scalar dPQ_dX1 = clientStrand == P ? s : -t;

        const Vec3x PQ = ( ps - qt);
        const Vec3x dDepl_dX0 = -( dPQ_dX0 * grad + dGrad_dX0.transpose() * PQ );
        const Vec3x dDepl_dX1 = -( dPQ_dX1 * grad + dGrad_dX1.transpose() * PQ );

        Mat3x dd_gradPQ_dX_dX;
        dd_gradPQ_dX_dX.setZero();
        for ( unsigned i = 0; i < 3; ++i )
        {
            dd_gradPQ_dX_dX += dGrad_dX_dXi[i] * PQ[i];
        }

        Eigen::Matrix<Scalar, 6, 6> j = Eigen::Matrix<Scalar, 6, 6>::Zero();

        // Let's consider that  d dcross / dP = Id
        j.block < 3, 3 > ( 0, 0 ) = dDepl_dX0 * dDepl_dX0.transpose()
                + depl * ( dd_gradPQ_dX_dX + dPQ_dX0 * dGrad_dX0 + dPQ_dX0 * dGrad_dX0.transpose() );
        j.block < 3, 3 > ( 0, 3 ) =
                dDepl_dX0 * dDepl_dX1.transpose()
                        + depl
                                * ( -dd_gradPQ_dX_dX + dPQ_dX0 * dGrad_dX1
                                        + dPQ_dX1 * dGrad_dX0.transpose() );
        j.block < 3, 3 > ( 3, 0 ) =
                dDepl_dX1 * dDepl_dX0.transpose()
                        + depl
                                * ( -dd_gradPQ_dX_dX + dPQ_dX1 * dGrad_dX0
                                        + dPQ_dX0 * dGrad_dX1.transpose() );
        j.block < 3, 3 > ( 3, 3 ) = dDepl_dX1 * dDepl_dX1.transpose()
                + depl * ( dd_gradPQ_dX_dX + dPQ_dX1 * dGrad_dX1 + dPQ_dX1 * dGrad_dX1.transpose() );

        //std::cout << clientStrand << " " << ( j - j.transpose() ).squaredNorm() << std::endl ;

        const unsigned vertexId = ( clientStrand == P ) ? m_edgeP : m_edgeQ;
        geometry.m_totalJacobian->edgeStencilAdd<6>( 4 * vertexId, -m_stiffness * j );

    }

#else
    const Scalar volume = tetraVolume();
    if ( volume * m_initialSign - m_thickness >= 0 )
    return;

    const Scalar correctedVolume = volume - m_initialSign * m_thickness;
    Eigen::Matrix<Scalar, 6, 6> hess;
    hess.block<3, 3>( 0, 0 ).setZero();
    hess.block<3, 3>( 3, 3 ).setZero();
    Eigen::Matrix<Scalar, 6, 1> grad;

    switch ( clientStrand )
    {
        case P:
        {
            hess.block<3, 3>( 0, 3 ) = hessP0P1TetraVolume();
            hess.block<3, 3>( 3, 0 ) = hess.block<3, 3>( 0, 3 ).transpose();

            grad.segment<3>( 0 ) = gradP0TetraVolume();
            grad.segment<3>( 3 ) = gradP1TetraVolume();

            Eigen::Matrix<Scalar, 6, 6> localJ = -2.0
            * ( hess * correctedVolume + grad * grad.transpose() ) * m_stiffness;

            geometry.m_totalJacobian->edgeStencilAdd<6>( 4 * m_edgeP, localJ );
        }
        break;

        case Q:
        {
            hess.block<3, 3>( 0, 3 ) = hessQ0Q1TetraVolume();
            hess.block<3, 3>( 3, 0 ) = hess.block<3, 3>( 0, 3 ).transpose();

            grad.segment<3>( 0 ) = gradQ0TetraVolume();
            grad.segment<3>( 3 ) = gradQ1TetraVolume();

            Eigen::Matrix<Scalar, 6, 6> localJ = -2.0
            * ( hess * correctedVolume + grad * grad.transpose() ) * m_stiffness;

            geometry.m_totalJacobian->edgeStencilAdd<6>( 4 * m_edgeQ, localJ );
        }
        break;

    }
#endif
}

void StrandStrandForce::accumulateCurrentF( VecXx& globalF, const ElasticStrand& strand )
{

    ClientNumber clientStrand;
    if ( &strand == m_strandP )
        clientStrand = P;
    else if ( &strand == m_strandQ )
        clientStrand = Q;
    else
        return; // Nothing to do if we are trying to apply this force to another strand. This should never happen...

    findTheVertices( strand.getCurrentState() );

    Scalar E;
    accumulateEF( clientStrand, globalF, E );

}

void StrandStrandForce::accumulateFutureE( Scalar& energy, const ElasticStrand& strand )
{
    ClientNumber clientStrand;
    if ( &strand == m_strandP )
        clientStrand = P;
    else if ( &strand == m_strandQ )
        clientStrand = Q;
    else
        return; // Nothing to do if we are trying to apply this force to another strand. This should never happen...

    findTheVertices( strand.getFutureState() );

    Scalar E;
    accumulateE( clientStrand, E );

}

void StrandStrandForce::accumulateE( ClientNumber clientStrand, Scalar& globalE ) const
{
#ifdef USE_DISTANCE_FOR_PENALTY_FORCES
    Scalar s, t;
    GetClosestPoints( s, t );

    if ( s == 0. || s == 1. || t == 0. || t == 1. )
        return;

    const Vec3x &ps = ( 1. - s ) * m_p0 + s * m_p1;
    const Vec3x &qt = ( 1. - t ) * m_q0 + t * m_q1;

    Vec3x edge, otherEdge;
    if ( clientStrand == P )
    {
        edge = ( m_p1 - m_p0 );
        otherEdge = ( m_q1 - m_q0 ).normalized();
    }
    else
    {
        edge = ( m_q1 - m_q0 );
        otherEdge = ( m_p1 - m_p0 ).normalized();
    }
    const Scalar edgeLen = edge.norm();

    Mat3x cross = -m_initialSign * crossMat( otherEdge );
    if ( clientStrand == Q )
    {
        cross *= -1;
    }

    const Vec3x grad = ( cross * edge ).normalized();
    const double rad = ( ps - qt ).dot( grad );
    const double depl = m_thickness - rad;

    if ( depl > 0. )
    {

        globalE += .5 * m_stiffness * std::pow( depl, 2 );

    }

#else
    const Scalar volume = tetraVolume();
    if ( volume * m_initialSign - m_thickness >= 0 )
    return;

    const Scalar correctedVolume = volume - m_initialSign * m_thickness;
    const Scalar localE = correctedVolume * correctedVolume * m_stiffness;

    globalE += localE;

#endif

}

void StrandStrandForce::accumulateEF( ClientNumber clientStrand, VecXx& globalF,
        Scalar& globalE ) const
{

#ifdef USE_DISTANCE_FOR_PENALTY_FORCES
    Scalar s, t;
    GetClosestPoints( s, t );

    if ( s == 0. || s == 1. || t == 0. || t == 1. )
        return;

    const Vec3x &ps = ( 1. - s ) * m_p0 + s * m_p1;
    const Vec3x &qt = ( 1. - t ) * m_q0 + t * m_q1;

    Vec3x edge, otherEdge;
    if ( clientStrand == P )
    {
        edge = ( m_p1 - m_p0 );
        otherEdge = ( m_q1 - m_q0 ).normalized();
    }
    else
    {
        edge = ( m_q1 - m_q0 );
        otherEdge = ( m_p1 - m_p0 ).normalized();
    }
    const Scalar edgeLen = edge.norm();

    Mat3x cross = -m_initialSign * crossMat( otherEdge );
    if ( clientStrand == Q )
    {
        cross *= -1;
    }

    const Vec3x grad = ( cross * edge ).normalized();
    const double rad = ( ps - qt ).dot( grad );
    const double depl = m_thickness - rad;

    if ( depl > 0. )
    {

        globalE += .5 * m_stiffness * std::pow( depl, 2 );

        const Mat3x dGrad_dX = ( Mat3x::Identity() - grad * grad.transpose() ) * cross / edgeLen;
        const Mat3x dGrad_dX0 = -dGrad_dX;
        const Mat3x dGrad_dX1 = dGrad_dX;

        // X = P|Q
        // alpha = d(Ps - Qt) / dX
        Scalar dPQ_dX0 = clientStrand == P ? ( 1 - s) : ( t - 1);
        Scalar dPQ_dX1 = clientStrand == P ? s : -t;

        const Vec3x dDepl_dX0 = -( dPQ_dX0 * grad + dGrad_dX0.transpose() * ( ps - qt ) );
        const Vec3x dDepl_dX1 = -( dPQ_dX1 * grad + dGrad_dX1.transpose() * ( ps - qt ) );

        const unsigned vertexId = ( clientStrand == P ) ? m_edgeP : m_edgeQ;
        globalF.segment < 3 > ( ( vertexId + 0 ) * 4 ) -= m_stiffness * depl * dDepl_dX0;
        globalF.segment < 3 > ( ( vertexId + 1 ) * 4 ) -= m_stiffness * depl * dDepl_dX1;

    }

#else
    const Scalar volume = tetraVolume();
    if ( volume * m_initialSign - m_thickness >= 0 )
    return;

    const Scalar correctedVolume = volume - m_initialSign * m_thickness;
    const Scalar localE = correctedVolume * correctedVolume * m_stiffness;

    globalE += localE;

    Eigen::Matrix<Scalar, 6, 1> localF;

    switch ( clientStrand )
    {
        case P:
        localF.segment<3>( 0 ) = -2.0 * gradP0TetraVolume() * correctedVolume * m_stiffness;
        localF.segment<3>( 3 ) = -2.0 * gradP1TetraVolume() * correctedVolume * m_stiffness;
        globalF.segment<3>( ( m_edgeP + 0 ) * 4 ) += localF.segment<3>( 0 );
        globalF.segment<3>( ( m_edgeP + 1 ) * 4 ) += localF.segment<3>( 3 );
        break;

        case Q:
        localF.segment<3>( 0 ) = -2.0 * gradQ0TetraVolume() * correctedVolume * m_stiffness;
        localF.segment<3>( 3 ) = -2.0 * gradQ1TetraVolume() * correctedVolume * m_stiffness;
        globalF.segment<3>( ( m_edgeQ + 0 ) * 4 ) += localF.segment<3>( 0 );
        globalF.segment<3>( ( m_edgeQ + 1 ) * 4 ) += localF.segment<3>( 3 );
        break;
    }

//    std::cout << "Client: " << clientStrand << " energy: " << localE << " force: " << localF
//            << '\n';
#endif

}

bool StrandStrandForce::isActive() const
{
//    return true ;
    return m_colliding;

    /*if ( ! SegmentsProjOnEachOther() ) return false ;
     //if ( ( tetraVolume() * m_initialSign - 3 * m_thickness ) > 0 ) return false ;

     Scalar s, t ;
     GetClosestPoints( s, t );
     const Vec3x &ps = ( 1. - s ) * m_p0 + s * m_p1;
     const Vec3x &qt = ( 1. - t ) * m_q0 + t * m_q1;

     return  ( ps - qt ).norm() < 4 * m_thickness  ;*/

}

void StrandStrandForce::updateCollisionData()
{
    //Collision in computed before both strands are resolved, in order to ensure symmetry

    VecXx globalF = VecXx::Zero( m_strandP->getNumVertices() * 4 );
    accumulateCurrentF( globalF, *m_strandP );

    Vec3x contactForce;
    contactForce = globalF.segment < 3 > ( ( m_edgeP + 0 ) * 4 ) + globalF.segment < 3
            > ( ( m_edgeP + 1 ) * 4 );

    if ( !isSmall( contactForce.squaredNorm() ) )
    {

        Scalar s, t;
        GetClosestPoints( s, t );
        const Vec3x &ps = ( 1. - s ) * m_p0 + s * m_p1;
        const Vec3x &qt = ( 1. - t ) * m_q0 + t * m_q1;
        const Scalar rad = ( ps - qt ).norm();

        StaticCollision &col = m_collision;

        col.force = contactForce;

        col.adhProps[StaticCollision::TANGENTIAL].coefficient = m_friction;
        col.adhProps[StaticCollision::TANGENTIAL].thickness = 10 * m_thickness;
        col.adhProps[StaticCollision::TANGENTIAL].perpThickness = m_thickness;

        col.adhProps[StaticCollision::NORMAL].coefficient =  m_stiffness
                * m_thickness;
        col.adhProps[StaticCollision::NORMAL].thickness = m_thickness;
        col.adhProps[StaticCollision::NORMAL].perpThickness = m_thickness;

        col.collisionPointData.resize( 15 );
        col.collisionPointData.segment < 3 > ( 0 ) = Vec3x( s, t, rad );

        col.collisionPointData.segment < 3 > ( 3 ) = m_p0;
        col.collisionPointData.segment < 3 > ( 6 ) = m_p1;
        col.collisionPointData.segment < 3 > ( 9 ) = m_q0;
        col.collisionPointData.segment < 3 > ( 12 ) = m_q1;

        col.localAbscissa = 0.;
        col.vertexId = 0;
        col.parentForce = this;

        if ( !m_colliding )
        {
            m_collisionAtCreation = m_collision;
            m_colliding = true;
        }

        //        if (clientStrand == P ) {
        //            std::cout << ps.transpose() << "  // " << getWorldCollisionPoint(col, strand, RuntimeCollisionForceBase::NORMAL).transpose() << std::endl ;
        //        } else {
        //            std::cout << qt.transpose() << "  // " << getWorldCollisionPoint(col, strand, RuntimeCollisionForceBase::NORMAL).transpose() << std::endl ;
        //        }
    }
    else
    {
        m_colliding = false;
        m_strandP_col_initialized = false;
        m_strandQ_col_initialized = false;
    }

}

void StrandStrandForce::registerCollision( const ElasticStrand& strand,
        StrandState& geometry ) const
{
    if( !m_colliding  ) return ;

    ClientNumber clientStrand;
    if ( &strand == m_strandP )
        clientStrand = P;
    else if ( &strand == m_strandQ )
        clientStrand = Q;
    else
        return; // Nothing to do if we are trying to apply this force to another strand. This should never happen...

    StaticCollision col;

    if ( &strand == m_strandP && !m_strandP_col_initialized )
    {
        col = m_collisionAtCreation;
        m_strandP_col_initialized = true;
    }
    else if ( &strand == m_strandQ && !m_strandQ_col_initialized )
    {
        col = m_collisionAtCreation;
        m_strandQ_col_initialized = true;
    }
    else
    {
        col = m_collision;
    }

    if ( &strand == m_strandP )
    {
        col.localAbscissa = m_collisionAtCreation.collisionPointData( 0 );
        col.vertexId = m_edgeP;

    }
    else
    {
        col.force *= -1;
        col.localAbscissa = m_collisionAtCreation.collisionPointData( 1 );
        col.vertexId = m_edgeQ;

    }

    geometry.staticCollisions().registerCollision( col );
}

Vec3x StrandStrandForce::getWorldCollisionNormal( const StaticCollision &collision,
        const ElasticStrand& strand ) const
{
    //findTheVertices( strand.getCurrentState() );

    return RuntimeCollisionForceBase::getWorldCollisionNormal( collision, strand );
}

Vec3x StrandStrandForce::getWorldCollisionPoint( const StaticCollision &collision,
        const ElasticStrand& strand ) const
{
    //findTheVertices( strand.getCurrentState() );

    Scalar s = collision.collisionPointData( 0 );
    Scalar t = collision.collisionPointData( 1 );

//    const Vec3x &p0 = m_collision.collisionPointData.segment< 3 >(3) ;
//    const Vec3x &p1 = m_collision.collisionPointData.segment< 3 >(6) ;
//    const Vec3x &q0 = m_collision.collisionPointData.segment< 3 >(9) ;
//    const Vec3x &q1 = m_collision.collisionPointData.segment< 3 >(12) ;
    const Vec3x& p0 = m_strandP->getVertex( m_edgeP );
    const Vec3x& p1 = m_strandP->getVertex( m_edgeP + 1 );
    const Vec3x& q0 = m_strandQ->getVertex( m_edgeQ );
    const Vec3x& q1 = m_strandQ->getVertex( m_edgeQ + 1 );

    Scalar rad = 0.;
    if ( collision.subspace == StaticCollision::NORMAL )
    {
        //s = m_collision.collisionPointData( 0 ) ;
        //t = m_collision.collisionPointData( 1 ) ;
        rad = m_collisionAtCreation.collisionPointData( 2 );
    }

    const Vec3x &ps = ( 1. - s ) * p0 + s * p1;
    const Vec3x &qt = ( 1. - t ) * q0 + t * q1;

    return ( ( &strand == m_strandP ) ? qt : ps ) + rad * collision.force.normalized();

}

bool StrandStrandForce::SegmentsProjOnEachOther() const
{
    const Vec3x& p0 = m_strandP->getVertex( m_edgeP );
    const Vec3x& p1 = m_strandP->getVertex( m_edgeP + 1 );
    const Vec3x& q0 = m_strandQ->getVertex( m_edgeQ );
    const Vec3x& q1 = m_strandQ->getVertex( m_edgeQ + 1 );

    Scalar s, t; //dummy
    return strandsim::SegmentsProjOnEachOther( p0, p1, q0, q1, s, t );
}

void StrandStrandForce::GetClosestPoints( Scalar &s, Scalar &t ) const
{
    strandsim::SegmentsProjOnEachOther( m_p0, m_p1, m_q0, m_q1, s, t );

    s = std::min( ( Scalar ) 1., std::max( s, ( Scalar ) 0. ) );
    t = std::min( ( Scalar ) 1., std::max( t, ( Scalar ) 0. ) );
}

std::pair<int, int> StrandStrandForce::SegmentsSituation() const
{
    const Vec3x& p0 = m_strandP->getVertex( m_edgeP );
    const Vec3x& p1 = m_strandP->getVertex( m_edgeP + 1 );
    const Vec3x& q0 = m_strandQ->getVertex( m_edgeQ );
    const Vec3x& q1 = m_strandQ->getVertex( m_edgeQ + 1 );

    return strandsim::SegmentsSituation( p0, p1, q0, q1 );
}

Scalar StrandStrandForce::getInitialSign() const
{
    return m_initialSign;
}

int& StrandStrandForce::getEdgeP()
{
    return m_edgeP;
}

int& StrandStrandForce::getEdgeQ()
{
    return m_edgeQ;
}

const ElasticStrand *StrandStrandForce::getStrandP() const
{
    return m_strandP;
}

const ElasticStrand *StrandStrandForce::getStrandQ() const
{
    return m_strandQ;
}

const ElasticStrand *StrandStrandForce::getOther( const ElasticStrand * me ) const
{
    if ( m_strandQ == me )
        return m_strandP;
    if ( m_strandP == me )
        return m_strandQ;

    return NULL;
}
