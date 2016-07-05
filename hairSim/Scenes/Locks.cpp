#include "Locks.h"

#include <random>
#define PI 3.14159265
using namespace std;

Locks::Locks():
Scene("Locks", "L Locks of hair forming an N braid")
{
    AddOption("lock_radius", "", 1.6 );

    GetScalarOpt("density") = 1.32;
    GetScalarOpt("youngs-modulus") = 3.9e+09;
    GetScalarOpt("viscosity") = 5e+8;
    GetScalarOpt("stretchingThreshold") = 1.;
    GetScalarOpt("radius") = 0.0037;
    GetIntOpt("maxNewtonIterations") = 30;

    AddOption("end_time","", 10.015 );
    // GetBoolOpt("render") = false;
    GetIntOpt("numberOfThreads") = 4;

    // Pre-setup to default values:
    GetScalarOpt("stochasticPruningFraction" ) = 0.5;

    AddOption("hairtie", "include a hair tie", false);
    GetBoolOpt("useProxRodRodCollisions") = true;
    GetBoolOpt("useCTRodRodCollisions") = true;
    GetScalarOpt("collisionRadius") = 0.01;
    GetScalarOpt("strand_mu") = 0.3;
    GetScalarOpt("dt") = 1e-2;

    GetScalarOpt("mesh_mu") = 0.1;
    GetStringOpt("checkpointDir") = "/Users/henrique/Desktop/LearnHair/build/Apps/StrandSimulator";
}

Locks::~Locks()
{}

void Locks::layeredLockPacking( std::vector< std::vector< Vec3 > >& rodPos )
{
    // create packed disk/circle to be repeated for each lock
    const unsigned L = 4; // layers (including center) of rings of strands
    const double dThetaMin = 3; // increment theta in search for non-collision
    const double R = GetScalarOpt("collisionRadius");
    const double b = 0.01 * R; // buffer of space between rods

    std::cout << "collision thickness: " << GetScalarOpt("collisionRadius") << std::endl;

    std::vector< Vec3 > centroids;
    centroids.push_back( Vec3::Zero() );

    double deg, degOffset, r;
    for( unsigned l = 1; l < L; ++l )
    {
        r = l * ( 2*R + b );
        deg = 0.0;
        degOffset = dThetaMin;
        bool firstInLayer = true;
        while( deg < 360.0 )
        { // attempt to add strand

            Scalar x = r * cos( deg * M_PI/180.0 );
            Scalar z = r * sin( deg * M_PI/180.0 );
            Vec3 potentialCentroid( x, 0.0, z );

            bool noOverlaps = true;
            for( unsigned c = 0; c < centroids.size(); ++c )
            { // check for overlaps
                if( (centroids[c] - potentialCentroid).norm() <= (2*R + b) ){
                    noOverlaps = false;
                    deg += dThetaMin;
                    break;
                }
            }

            if( noOverlaps ){
                if( firstInLayer ){
                    firstInLayer = false;
                    degOffset = deg;
                }
                centroids.push_back( potentialCentroid );
                deg += degOffset;
            }
        }
    }

    unsigned rods_per_lock = centroids.size();
    double outermostR = centroids.back().norm();
    std::cout << "outermost radius of lock: " << outermostR << std::endl;
    std::cout << "number of strands per lock: " << rods_per_lock << std::endl;

    if( outermostR == 0.0 ) outermostR = GetScalarOpt("collisionRadius") + 0.01;

    const unsigned N = 3; // number of centerlines
    const unsigned M = 2; // even number < N
    const double A = 2.05 * outermostR; // scaling for X
    const double B = 2.05 * outermostR; // scaling for Z
    const double t = 1.0 * 0.5; // change in Y per edge
    const unsigned numVerts = 20; // length of strands
    rodPos.resize( rods_per_lock * N );

    Scalar y = 0.0;
    for( unsigned v = 0; v < numVerts; ++v )
    {
        for( unsigned rod = 0; rod < N; ++rod )
        {
            Scalar x = A * sin( y + rod * 2 * PI / N );
            Scalar z = B * sin( M * ( y + rod * 2 * PI / N ) );

            for( unsigned strand = 0; strand < rods_per_lock; ++strand )
            {
                Vec3 vertNext = Vec3( x, y, z );
                vertNext += centroids[strand];
                rodPos[rod * rods_per_lock + strand].push_back( vertNext );                
            }
        }
        y -= t;
    }    

}

int hairtieNum = -1;
bool includeLockTie = false;
void Locks::includeHairTie()
{
    includeLockTie = true;

    // rod options
    Scalar radiusA = GetScalarOpt("radius");
    Scalar youngsModulus = GetScalarOpt("youngs-modulus");
    Scalar shearModulus = GetScalarOpt("shear-modulus");
    Scalar density = GetScalarOpt("density");
    Scalar airDrag = GetScalarOpt("air-drag");
    Scalar viscosity = GetScalarOpt("viscosity");
    Scalar baseRotation = GetScalarOpt("base-rotation");

    double depth = -7.5;
    double radius = 0.2;
    vector< Vec3 > vertices;
    vertices.push_back( Vec3( radius, depth, 0.0 ) );
    vertices.push_back( Vec3( 0.0, depth, radius ) );
    vertices.push_back( Vec3( -radius, depth, 0.0 ) );
    vertices.push_back( Vec3( 0.0, depth, -radius ) );
    vertices.push_back( Vec3( radius, depth, 0.0 ) );

    int num_DoFs = 4 * vertices.size() - 1;
    VecXd dofs( num_DoFs );
    for ( unsigned i = 0; i < vertices.size(); ++i ){
        dofs.segment<3>( i * 4 ) = vertices[i];
    }

    DOFScriptingController* controller = new DOFScriptingController( );

    ElasticStrandParameters* params = new ElasticStrandParameters( radiusA, youngsModulus, shearModulus, density, viscosity, airDrag, baseRotation );
    ElasticStrand* strand = new ElasticStrand( dofs, *params, controller, 0.15 );
    strand->setGlobalIndex( m_strands.size() );
    hairtieNum = m_strands.size();
    setRodCollisionParameters( *strand );

    m_strands.push_back( strand );

    Vec3 zero( 0.,0.,0. );
    ElasticStrand* hairband = m_strands[ m_strands.size() - 1 ];
    Mat3x identity = Mat3x::Identity();
    double push = 0.5;
    Vec3 translate = Vec3( push, 0.0, 0.0 );
    SceneUtils::transformRodRootVtx( hairband, identity, zero, translate , 0 );
    translate = Vec3( 0.0, 0.0, push );
    SceneUtils::transformRodRootVtx( hairband, identity, zero, translate, 1 );
    translate = Vec3( -push, 0.0, 0.0 );
    SceneUtils::transformRodRootVtx( hairband, identity, zero, translate, 2 );
    translate = Vec3( 0.0, 0.0, -push );
    SceneUtils::transformRodRootVtx( hairband, identity, zero, translate, 3 );
    translate = Vec3( push, 0.0, 0.0 );
    SceneUtils::transformRodRootVtx( hairband, identity, zero, translate, 4 );
}

void Locks::setupStrands()
{
    // rod options
    Scalar radiusA = GetScalarOpt("radius");
    Scalar youngsModulus = GetScalarOpt("youngs-modulus");
    Scalar shearModulus = GetScalarOpt("shear-modulus");
    Scalar density = GetScalarOpt("density");
    Scalar airDrag = GetScalarOpt("air-drag");
    Scalar viscosity = GetScalarOpt("viscosity");
    Scalar baseRotation = GetScalarOpt("base-rotation");

    // sample rod positions
    vector< vector< Vec3 > >strands;
    layeredLockPacking( strands );    

    // Vec3 CoM = Vec3(0.0,0.0,0.0);
    int nVertices;

    std::cout << "numStrands : " << strands.size() << std::endl;
    // more rod params
    int num_DoFs;
    int rod_id = 0;
    for( unsigned i = 0; i < strands.size(); ++i )
    {
        vector< Vec3 > vertices = strands[i];
        num_DoFs = 4 * vertices.size() - 1;
        nVertices = vertices.size();

        // initial strand position
        VecXd dofs( num_DoFs );
        for ( unsigned i = 0; i < vertices.size(); ++i ){
            dofs.segment<3>( i * 4 ) = vertices[i];
        }

        DOFScriptingController* controller = new DOFScriptingController( );
        controller->freezeVertices( 0, true );

        ElasticStrandParameters* params = new ElasticStrandParameters( 
                radiusA, youngsModulus, shearModulus, 
                density, viscosity, airDrag, baseRotation );
        ElasticStrand* strand = new ElasticStrand( dofs, *params, controller, GetScalarOpt("collisionRadius") );
        strand->setGlobalIndex( rod_id );
        setRodCollisionParameters( *strand );

        Vec2Array kappas = strand->alterRestKappas();
        for( unsigned i = 0; i < kappas.size(); ++i )
        {   // make strands straight regardless of starting configuration
            strand->alterRestKappas()[i].setZero();
        }

        m_strands.push_back( strand );
        ++rod_id;
    }

    if( GetBoolOpt("hairtie") ) includeHairTie();

    std::cout << "# num strands: \n" << m_strands.size() <<'\n';
    std::cout << "# num dofs per strand: \n" << 4 * strands[0].size() - 1 <<'\n';
}

void Locks::setupMeshes()
{}

bool Locks::executeScript()
{
    if ( getTime() >= GetScalarOpt("end_time") )
    {
        std::cout << "# Simulation complete. Exiting." << std::endl;
        exit(0);
    }

    Vec3 zero( 0.,0.,0. );
    if( getTime() == 0.0 + m_dt && includeLockTie )
    {
        std::cout << "expanding hairtie" << std::endl;
        ElasticStrand* hairband = m_strands[ m_strands.size() - 1 ];
        Mat3x identity = Mat3x::Identity();
        double push = -0.0;
        Vec3 translate = Vec3( push, 0.0, 0.0 );
        SceneUtils::transformRodRootVtx( hairband, identity, zero, translate , 0 );
        translate = Vec3( 0.0, 0.0, push );
        SceneUtils::transformRodRootVtx( hairband, identity, zero, translate, 1 );
        translate = Vec3( -push, 0.0, 0.0 );
        SceneUtils::transformRodRootVtx( hairband, identity, zero, translate, 2 );
        translate = Vec3( 0.0, 0.0, -push );
        SceneUtils::transformRodRootVtx( hairband, identity, zero, translate, 3 );
        translate = Vec3( push, 0.0, 0.0 );
        SceneUtils::transformRodRootVtx( hairband, identity, zero, translate, 4 );
    }

    if( getTime() == 0.0 + 2 * m_dt && includeLockTie )
    {
        ElasticStrand* hairband = m_strands[ m_strands.size() - 1 ];
        hairband->dynamics().getScriptingController()->clear();
    }
    
    // // freezeTriangleObject( *currentMesh );
    // unsigned rod = 0;
    // for(auto rd_itr = m_strands.begin(); rd_itr != m_strands.end(); ++rd_itr, ++rod )
    // {
    //     if( rod == m_strands.size() - 1 ) continue;

    //     freezeVertex( **rd_itr, 0);
    //     // if( getTime() >= 0.25 ){
    //     //     freezeVertex( **rd_itr, nVertices - 1);
    //     // }
    // }

    bool script = false;
    int scriptType = 2;
    if ( script)
    {
        if( scriptType == 1 )
        {
        }
        else if( scriptType == 2 )
        {
            double swayRate = 20 * M_PI;

            Mat3x identity = Mat3x::Identity();
            // std::cout << "script cos: " << cos( swayRate * getTime() ) << std::endl;
            Vec3 translate = Vec3( 0.1 * sin( swayRate * getTime() ), 0.0, 0.1 * cos( swayRate * getTime() ) );

            int rodCount = 0;
            for(auto rd_itr = m_strands.begin(); rd_itr != m_strands.end(); ++rd_itr, ++rodCount)
            {
                if( rodCount == hairtieNum ) continue; // hairtie
                SceneUtils::transformRodRootVtx( *rd_itr, identity, zero, translate, 0 );
            }
        }
        else if( scriptType == 3 && getTime() > 0.0 + 5000 * m_dt  && includeLockTie )
        { // Move hairtie
            double swayRate = 20 * M_PI;

            Mat3x identity = Mat3x::Identity();
            // std::cout << "script cos: " << cos( swayRate * getTime() ) << std::endl;
            Vec3 translate = Vec3( 0.1 * sin( swayRate * getTime() ), 0.001, 0.1 * cos( swayRate * getTime() ) );

            int rodCount = 0;
            for(auto rd_itr = m_strands.begin(); rd_itr != m_strands.end(); ++rd_itr, ++rodCount)
            {
                if( rodCount == 20 )
                {
                    for( unsigned i = 0; i < (*rd_itr)->getNumVertices(); ++i ){
                        SceneUtils::transformRodRootVtx( *rd_itr, identity, zero, translate, i );
                    }
                }
            }
        }

    }
    else{
      // freezeTriangleObject( *currentMesh );
      unsigned rod = 0;
      for(auto rd_itr = m_strands.begin(); rd_itr != m_strands.end(); ++rd_itr, ++rod )
      {
        if( rod == m_strands.size() - 1 ) continue;

        // freezeRodRoot( **rd_itr );
        SceneUtils::freezeVertex( *rd_itr, 0);
        // if( getTime() >= 0.25 ){
        //     freezeVertex( **rd_itr, nVertices - 1);
        // }
      }
    }


    return true;  
}

