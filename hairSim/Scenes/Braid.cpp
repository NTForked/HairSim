#include "Braid.hh"
#include <random>

#define PI 3.14159265
using namespace std;

Braid::Braid() :
ProblemStepper("Braid", "N Locks of hair forming a braid ")
{

    AddOption("num_nurbs", "", 1500 );
    AddOption( "directory" , "file with braid maya nurbs file" , "/Users/henrique/Desktop/DutchBraid" );

    AddOption("totalLength","enforced strand length", 30.0 );
    AddOption("sphere_mesh_filename","","assets/TriangulatedSphere.obj"); 
    AddOption("geometryscale","uniform rescaling of sphere and hair root locations", 0.3 );


    AddOption("scripting_on","", true );
    AddOption("rotatescale","uniform rescaling of rotational speeds", 1.0);
    AddOption("time_offset","start scripting at time t + offset", 1.75);


    AddOption("lock_diameter", "", 0.12 ); // lock diameter and swap_proximity must be tuned together to look right...
    AddOption("swap_proximity", "", 0.5 );
    AddOption("rods_per_lock", "", 1 ); // rods per lock as well
    AddOption("swap_partitions", "", 5 );
    AddOption("swap_count", "", 7 );
    AddOption("y_offset", "", -1.0 );

    GetScalarOpt("density") = 1.32;
    GetScalarOpt("youngs-modulus") = 3.9e+09;
    GetScalarOpt("viscosity") = 5e+8;
    GetScalarOpt("stretchingThreshold") = 1.;
    GetScalarOpt("radiusA") = 0.0037;
    GetScalarOpt("radiusB") = 0.0037;
    GetIntOpt("maxNewtonIterations") = 30;

    AddOption("end_time","", 20.015 );  

    // GetBoolOpt("render") = false;
    GetIntOpt("numberOfThreads") = 1;

    // Pre-setup to default values:
    GetScalarOpt( "stochasticPruningFraction" ) = 0.5;
    GetBoolOpt("useProxRodRodCollisions") = true;
    GetBoolOpt("useCTRodRodCollisions") = true;
    GetScalarOpt("percentCTRodRodCollisionsAccept") = 100.0;
    GetScalarOpt("selfCollisionsRadius") = 0.16; 


    GetScalarOpt("dt") = 1e-2;
    GetScalarOpt("strand_mu") = 0.7;
    GetScalarOpt("mesh_mu") = 0.4;
    GetStringOpt("checkpointDir") = "/Users/henrique/Desktop/LearnHair/build/Apps/StrandSimulator";

}

Braid::~Braid()
{    
}

void Braid::generateSinusoidalBraidVertices( vector< vector< Vec3x > >& rodPos )
{
    // This formula for odd-locked braids 
    // from [Capturing Braided Hairstyles - Hu et al. 2014]

    const unsigned N = 5; // number of centerlines
    const unsigned M = 2; // even number < N
    const double A = 1.0 * 0.333; // scaling for X
    const double B = 1.0 * 0.667; // scaling for Z
    const double t = 1.0 * 0.5; // change in Y per edge
    const unsigned numVerts = 20; // length of strands

    rodPos.resize( N );
    Scalar y = 0.0;
    for( unsigned v = 0; v < numVerts; ++v )
    {
        for( unsigned rod = 0; rod < rodPos.size(); ++rod )
        {
            Scalar x = A * sin( y + rod * 2 * PI / N );
            Scalar z = B * sin( M * ( y + rod * 2 * PI / N ) );
            // std::cout << "rod: " << rod << " x: " << x << " y: " << y << " z: " << z << std::endl;
            Vec3x vertNext = Vec3x( x, y, z );
            rodPos[rod].push_back( vertNext );
        }
        y -= t;
    }
}

void Braid::generateBraidVertices( vector< vector< Vec3x > >& rodPos )
{
    //// Currently always has 3 locks generated on the XZ plane, aligned with X and hanging down in Y

    // First generate starting points for each lock
    double lockRadius = GetScalarOpt("lock_diameter") / 2.0;
    std::uniform_real_distribution<double> radius( 0.0, lockRadius ); // [min, max)
    std::uniform_real_distribution<double> degree( 0.0, 360.0 );  // [min, max)

    std::mt19937 rng; // Mersenne Twister: Good quality random number generator
    rng.seed(); //Initialize with default, (deterministic?) seed

    // Scalar y = GetScalarOpt("y_offset");
    Scalar y = - GetScalarOpt("geometryscale") * 4.0; //cm


    unsigned numRodsPerLock = GetIntOpt("rods_per_lock");
    // vector< vector <Vec3x> > 
    rodPos.resize( 3 * numRodsPerLock );
    for( unsigned lock = 0; lock < 3; ++lock ){
        for( unsigned rodInLock = 0; rodInLock < numRodsPerLock; ++rodInLock ){
            // Generate r and theta value (polar coords)
            Scalar rad = radius(rng);
            Scalar theta = degree(rng);

            // Convert to Cartesian
            Scalar x = rad * cos( theta * M_PI/180.0 );
            Scalar z = rad * sin( theta * M_PI/180.0 );
            x += (-2.0 * lockRadius) + lock * (2*lockRadius); // shifting to center of lock

            // Add to Rods list, grouped by lock
            rodPos[ lock * numRodsPerLock + rodInLock ].push_back( Vec3x( x, y, z ) );
        }
    }

    Scalar dY =  - GetScalarOpt("swap_proximity") / GetIntOpt("swap_partitions");

    for( unsigned rod = 0; rod < rodPos.size(); ++rod ){
        Vec3x vertNext = rodPos[rod].back();

        vertNext += Vec3x( 0.0, dY, 0.0 );
        rodPos[rod].push_back( vertNext );
    }


    Scalar dTheta = ( 360.0 / GetIntOpt("swap_partitions") ) * M_PI / 180.0;
    Mat3x rotY;
    rotY <<
    cos( dTheta ),  0, sin( dTheta ),
        0,          1,      0,
    -sin( dTheta ), 0, cos( dTheta );

    bool swapLeft = true; // Alternate between swapping left and right with middle
    Scalar tangentShift;
    Vec3x pivot;
    unsigned startIndex;
    unsigned endIndex;
    // Number of times an outside lock is swapped with the middle lock
    for( int swapCount = 0; swapCount < GetIntOpt("swap_count"); ++swapCount ){

        if( swapLeft ){
            tangentShift = -lockRadius;
            startIndex = 0;
            endIndex = 2 * numRodsPerLock;
        }
        else{
            tangentShift = lockRadius;
            startIndex = numRodsPerLock;
            endIndex = 3 * numRodsPerLock;
        }

        // divide this swap into an appropriate number of swap partitions for smoothness
        for( int rotationCount = 0; rotationCount < GetIntOpt("swap_partitions"); ++rotationCount ){

            pivot = Vec3x( tangentShift, y, 0.0 ); // rotate about previous level
            // rotate locks of interest (and their rod points) about tangent point of lock circles
            for( unsigned rod = 0; rod < rodPos.size(); ++rod ){
                Vec3x vertNext = rodPos[rod].back();

                if( startIndex <= rod && rod < endIndex ){
                    vertNext = rotY * (vertNext - pivot) + pivot;
                }
                vertNext += Vec3x( 0.0, dY, 0.0 );
                rodPos[rod].push_back( vertNext );
            }
            y += dY; // but keep strands moving down
        }
        swapLeft = !swapLeft;   
    }
}

bool includeTie = false;

void Braid::includeHairTie()
{
    includeTie = true;

/*
    GetScalarOpt("density") = 1.32;
    GetScalarOpt("youngs-modulus") = 3.9e+09;
    GetScalarOpt("viscosity") = 5e+8;
    GetScalarOpt("stretchingThreshold") = 1.;
    GetScalarOpt("selfCollisionsRadius") = 0.0025;
    GetScalarOpt("radiusA") = 0.0037;
    GetScalarOpt("radiusB") = 0.0037;
*/

    // rod options
    Scalar radiusA = 0.01; //GetScalarOpt("radiusA");
    Scalar radiusB = 0.01; //GetScalarOpt("radiusB");
    Scalar youngsModulus = GetScalarOpt("youngs-modulus");
    Scalar shearModulus = GetScalarOpt("shear-modulus");
    Scalar density = GetScalarOpt("density");
    Scalar airDrag = GetScalarOpt("air-drag"); //1e-4
    Scalar viscosity = GetScalarOpt("viscosity");
    Scalar baseRotation = GetScalarOpt("base-rotation");

    double strandRadius = 0.15;
    double depth = -7.5;
    double radius = 1.1;
    double dOffset = 0.0; // 0.25
    vector< Vec3x > vertices;
    vertices.push_back( Vec3x( radius, depth, 0.0 ) );
    vertices.push_back( Vec3x( 0.0, depth, radius ) );
    vertices.push_back( Vec3x( -radius, depth, 0.0 ) );
    vertices.push_back( Vec3x( 0.0, depth, -radius ) );
    vertices.push_back( Vec3x( radius, depth + dOffset, 0.0 ) );
    // vertices.push_back( Vec3x( 0.0, depth + dOffset, radius ) );

    int num_DoFs = 4 * vertices.size() - 1;
    VecXd dofs( num_DoFs );
    for ( unsigned i = 0; i < vertices.size(); ++i ){
        dofs.segment<3>( i * 4 ) = vertices[i];
    }

    // CoM += vertices.back();
    strandsim::Vec3xArray scripted_vertices;
    scripted_vertices.push_back( vertices[0] );
    scripted_vertices.push_back( vertices[ vertices.size() - 1] );
    DOFScriptingController* controller = new DOFScriptingController( scripted_vertices );

    ElasticStrandParameters* params = new ElasticStrandParameters( radiusA, radiusB, youngsModulus, shearModulus, density, viscosity, airDrag, baseRotation );
    ElasticStrand* strand = new ElasticStrand( dofs, *params, controller, strandRadius );
    strand->setGlobalIndex( m_rodDatum.size() );
    setRodCollisionParameters( *strand );
    strand->collisionParameters().m_frictionCoefficient = 1.55;

    m_strands.push_back( strand );
    
    // extra stuff for render, etc...
    RodData* rd = new RodData( *strand, *controller );
    m_rodDatum.push_back( rd );

    Vec3x zero( 0.,0.,0. );
    RodData* hairband = m_rodDatum[ m_rodDatum.size() - 1 ];
    Mat3x identity = Mat3x::Identity();
    double push = 0.8;
    Vec3x translate = Vec3x( push, 0.0, 0.0 );
    transformRodRootVtx( *hairband, identity, zero, translate , 0 );
    translate = Vec3x( 0.0, 0.0, push );
    transformRodRootVtx( *hairband, identity, zero, translate, 1 );
    translate = Vec3x( -push, 0.0, 0.0 );
    transformRodRootVtx( *hairband, identity, zero, translate, 2 );
    translate = Vec3x( 0.0, 0.0, -push );
    transformRodRootVtx( *hairband, identity, zero, translate, 3 );
    translate = Vec3x( push, 0.0, 0.0 );
    transformRodRootVtx( *hairband, identity, zero, translate, 4 );
    // translate = Vec3x( 0.0, 0.0, push );
    // transformRodRootVtx( *hairband, identity, zero, translate, 5 );
}


int nVertices;
static bool nurbs = false;
void Braid::setupStrands()
{

    if( nurbs ){
        loadNurbs();
        return;
    }

    // rod options
    Scalar radiusA = GetScalarOpt("radiusA");
    Scalar radiusB = GetScalarOpt("radiusB");
    Scalar youngsModulus = GetScalarOpt("youngs-modulus");
    Scalar shearModulus = GetScalarOpt("shear-modulus");
    Scalar density = GetScalarOpt("density");
    Scalar airDrag = GetScalarOpt("air-drag");
    Scalar viscosity = GetScalarOpt("viscosity");
    Scalar baseRotation = GetScalarOpt("base-rotation");

    double strandRadius = 0.15;
    const Scalar totalLength = GetScalarOpt("totalLength");

    // sample rod positions
    vector< vector< Vec3x > >strands;
    // generateBraidVertices( strands ); 
    generateSinusoidalBraidVertices( strands );    

    // Vec3x CoM = Vec3x(0.0,0.0,0.0);

    // more rod params
    int num_DoFs;
    int rod_id = 0;
    for( unsigned i = 0; i < strands.size(); ++i )
    {
        vector< Vec3x > vertices = strands[i];
        num_DoFs = 4 * vertices.size() - 1;

        nVertices = vertices.size();

        Scalar length = 0.0;
        for ( unsigned i = 0; i < vertices.size() - 1; i++ )
            length += ( vertices[i + 1] - vertices[i] ).norm();
        std::cout << "actual_length: " << length <<  std::endl;
        // Enforce the total length
        // for ( unsigned i = 0; i < vertices.size(); i++ )
        //     vertices[i] *= totalLength / length;

        // initial strand position
        VecXd dofs( num_DoFs );
        for ( unsigned i = 0; i < vertices.size(); ++i ){
            dofs.segment<3>( i * 4 ) = vertices[i];
        }

        // CoM += vertices.back();
        strandsim::Vec3xArray scripted_vertices;
        scripted_vertices.push_back( vertices[0] );
        scripted_vertices.push_back( vertices[ nVertices - 1] );
        DOFScriptingController* controller = new DOFScriptingController( scripted_vertices );
        controller->freezeVertices( 0, true );
        // controller->freezeVertices( vertices.size() - 1 );


        ElasticStrandParameters* params = new ElasticStrandParameters( radiusA, radiusB, youngsModulus, shearModulus, density, viscosity, airDrag, baseRotation );
        ElasticStrand* strand = new ElasticStrand( dofs, *params, controller, strandRadius );
        strand->setGlobalIndex( rod_id );
        setRodCollisionParameters( *strand );

        Vec2xArray kappas = strand->alterRestKappas();
        for( unsigned i = 0; i < kappas.size(); ++i ){
            strand->alterRestKappas()[i].setZero();
        }

        m_strands.push_back( strand );
        
        // extra stuff for render, etc...
        RodData* rd = new RodData( *strand, *controller );
        m_rodDatum.push_back( rd );
        ++rod_id;
    }

    includeHairTie();

    std::cout << "# num strands: \n" << m_strands.size() <<'\n';
    std::cout << "# num dofs per strand: \n" << 4 * strands[0].size() - 1 <<'\n';

}

void Braid::setupMeshes()
{
    SimpleMeshController* mesh_controller = new SimpleMeshController( 0., m_dt );
    m_meshScripting_controllers.push_back( mesh_controller );
}

void Braid::loadNurbs()
{

    int numStrands = 10;
    // int numStrands = GetIntOpt("num_nurbs");
    std::string inputdirectory = GetStringOpt("directory");
    std::stringstream name;
    name << inputdirectory << "/dutch_" << GetIntOpt("num_nurbs") << ".ma";
    
    std::ifstream infile( name.str().c_str() );
    
    if( !infile.is_open() )
    {
        std::cerr << "\033[31;1mPlayback:\033[m can't load " << name << std::endl;
        exit(1);
    }
    
    //jump over maya nurbs header
    std::string header_string;
    for(int i=0; i < 1; ++i) getline( infile, header_string );
    
    bool first = true;

    int rod_id = 0;
    int numVerts;
    double x, y, z;
    while( !infile.eof() )
    {
        //skip header
        for( int i = 0; i < 7; ++i ){
            getline( infile, header_string );
        }
        
        //get number of verts:
        getline( infile, header_string );
        std::istringstream vertCount( header_string );
        vertCount >> numVerts;

        int delayedStart = numVerts / 2;

        std::vector<Vec3x> vertices;
        for( int v = 0; v < numVerts; ++v )
        {
            std::string vert_string;
            getline( infile, vert_string );

            // if( numVerts > 200 && v < 90 ) continue; // skip to the braid
            // if( numVerts > 150 && v < 60 ) continue; // skip to the braid
            if( v < delayedStart ) continue;

            std::istringstream vertstream( vert_string );
            vertstream >> x >> y >> z;
            Vec3x pos( x, -y, z );
            vertices.push_back( pos );
        }

        if( first ){
            first = false;
            continue;
        }

        // rod options
        Scalar radiusA = GetScalarOpt("radiusA");
        Scalar radiusB = GetScalarOpt("radiusB");
        Scalar youngsModulus = GetScalarOpt("youngs-modulus");
        Scalar shearModulus = GetScalarOpt("shear-modulus");
        Scalar density = GetScalarOpt("density");
        Scalar airDrag = GetScalarOpt("air-drag");
        Scalar viscosity = GetScalarOpt("viscosity");
        Scalar baseRotation = GetScalarOpt("base-rotation");
        
        // initial strand position
        int num_DoFs = 4 * vertices.size() - 1;
        VecXd dofs( num_DoFs );
        for ( int i = 0; i < dofs.size(); i += 4 )
            dofs.segment<3>( i ) = vertices[i / 4];
        
        strandsim::Vec3xArray scripted_vertices;
        scripted_vertices.push_back( vertices[0] );
        // scripted_vertices.push_back( vertices[1] );
        DOFScriptingController* controller = new DOFScriptingController( scripted_vertices );
        controller->freezeVertices(0);
        
        ElasticStrandParameters* params = new ElasticStrandParameters( radiusA, radiusB, youngsModulus, shearModulus, density, viscosity, airDrag, baseRotation );
        ElasticStrand* strand = new ElasticStrand( dofs, *params, controller );
        strand->setGlobalIndex( rod_id );
        setRodCollisionParameters( *strand );

        Vec2xArray kappas = strand->alterRestKappas();
        std::cout << "K_size " << kappas.size() <<std::endl;
        for( unsigned i = 0; i < kappas.size(); ++i ){
            strand->alterRestKappas()[i].setZero();
        }

        m_strands.push_back( strand );
        
        // extra stuff for render, etc...
        RodData* rd = new RodData( *strand, *controller );
        m_rodDatum.push_back( rd );
        std::cout<< "loaded rod " << rod_id << " with verts = "<< numVerts << std::endl;

        ++rod_id;
        if( rod_id >= numStrands ) break;
    }
    
    // std::cout << "\033[35;1mPlayBack message:\033[m loaded " << m_rodDatum.size() << " rods. from: " << name.str() << std::endl;
    std::cout << "\033[35;1mPlayBack message:\033[m loaded " << m_rodDatum.size() << " rods " << std::endl;
}

bool Braid::executeScript()
{
    if ( getTime() >= GetScalarOpt("end_time") )
    {
        std::cout << "# Simulation complete. Exiting." << std::endl;
        exit(0);
    }

    Vec3x zero( 0.,0.,0. );

    if( getTime() == 0.0 && includeTie )
    {
        std::cout << "expanding hairtie" << std::endl;
    }
    if( getTime() == 0.0 + m_dt && includeTie )
    {
        RodData* hairband = m_rodDatum[ m_rodDatum.size() - 1 ];
        Mat3x identity = Mat3x::Identity();
        double push = -0.5;
        Vec3x translate = Vec3x( push, 0.0, 0.0 );
        transformRodRootVtx( *hairband, identity, zero, translate , 0 );
        translate = Vec3x( 0.0, 0.0, push );
        transformRodRootVtx( *hairband, identity, zero, translate, 1 );
        translate = Vec3x( -push, 0.0, 0.0 );
        transformRodRootVtx( *hairband, identity, zero, translate, 2 );
        translate = Vec3x( 0.0, 0.0, -push );
        transformRodRootVtx( *hairband, identity, zero, translate, 3 );
        translate = Vec3x( push, 0.0, 0.0 );
        transformRodRootVtx( *hairband, identity, zero, translate, 4 );
        // translate = Vec3x( 0.0, 0.0, push );
        // transformRodRootVtx( *hairband, identity, zero, translate, 5 );
    }
    if( getTime() == 0.0 + 2 * m_dt && includeTie )
    {
        RodData* hairband = m_rodDatum[ m_rodDatum.size() - 1 ];
        hairband->getDofController().m_scriptedDegreesOfFreedom.clear();
    }

    if( getTime() == 0.0 + 15 * m_dt && includeTie )
    {
        GetBoolOpt("useProxRodRodCollisions") = true;
    }

    if( nurbs ) return true; // nurbs only have first vert scripted, below script expects both strand endpoints

    
    double thetaxrate = GetScalarOpt("rotatescale") * 2. * M_PI;
    double thetayrate = GetScalarOpt("rotatescale") * 2. * M_PI;
    double thetazrate = GetScalarOpt("rotatescale") * 2. * M_PI;
    
    Mat3x rotx;
    rotx << 1,0,0,
    0,cos(thetaxrate*getDt()),-sin(thetaxrate*getDt()),
    0,sin(thetaxrate*getDt()),cos(thetaxrate*getDt());
    Mat3x rotxT = rotx.transpose();
    
    Mat3x roty;
    roty << cos(thetayrate*getDt()),0,sin(thetayrate*getDt()),
    0,1,0,
    -sin(thetayrate*getDt()),0,cos(thetayrate*getDt());
    Mat3x rotyT = roty.transpose();
    
    Mat3x rotz;
    rotz << cos(thetazrate*getDt()),-sin(thetazrate*getDt()),0,
    sin(thetazrate*getDt()),cos(thetazrate*getDt()),0,
    0,0,1;
    Mat3x rotzT = rotz.transpose();
    
    TriangularMesh* currentMesh = m_meshScripting_controllers[0]->getCurrentMesh();
    
    int scriptType = 3;
    if ( GetBoolOpt("scripting_on") )
    {
        if( scriptType == 1 )
        {
        }
        else if( scriptType == 2 )
        {
            double swayRate = 20 * M_PI;

            Mat3x identity = Mat3x::Identity();
            // std::cout << "script cos: " << cos( swayRate * getTime() ) << std::endl;
            Vec3x translate = Vec3x( 0.1 * sin( swayRate * getTime() ), 0.0, 0.1 * cos( swayRate * getTime() ) );

            int rodCount = 0;
            for(auto rd_itr = m_rodDatum.begin(); rd_itr != m_rodDatum.end(); ++rd_itr, ++rodCount)
            {
                if( rodCount == 5 ) continue; // hairtie
                transformRodRootVtx( **rd_itr, identity, zero, translate, 0 );
            }
        }
        else if( scriptType == 3 && getTime() > 0.0 + 5 * m_dt )
        { // Move hairtie
            double swayRate = 20 * M_PI;
            double scale = 1e-1;

            Mat3x identity = Mat3x::Identity();
            // std::cout << "script cos: " << cos( swayRate * getTime() ) << std::endl;
            Vec3x translate = Vec3x( scale * sin( swayRate * getTime() ), 0.0, scale * cos( swayRate * getTime() ) );

            int rodCount = 0;
            for(auto rd_itr = m_rodDatum.begin(); rd_itr != m_rodDatum.end(); ++rd_itr, ++rodCount)
            {
                if( rodCount == 5 )
                {
                    for( unsigned i = 0; i < (*rd_itr)->m_strand.getNumVertices(); ++i ){
                        transformRodRootVtx( **rd_itr, identity, zero, translate, i );
                    }
                }
            }
        }

    }
    else{
      freezeTriangleObject( *currentMesh );
      unsigned rod = 0;
      for(auto rd_itr = m_rodDatum.begin(); rd_itr != m_rodDatum.end(); ++rd_itr, ++rod )
      {
        if( rod == m_rodDatum.size() - 1 ) continue;

        // freezeRodRoot( **rd_itr );
        freezeVertex( **rd_itr, 0);
        // if( getTime() >= 0.25 ){
        //     freezeVertex( **rd_itr, nVertices - 1);
        // }
      }
    }

  
    return true;
  
}

