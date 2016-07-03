#include "Aleka.h"

Aleka::Aleka() :
Scene("Aleka", "Strands dragging orthogonally across other fixed strands"),
m_radius(3.)
{
    AddOption( m_problemName, m_problemDesc, "" );
    
    // Global opts
    GetScalarOpt("dt") = 1e-3;
    GetVecOpt("gravity") = Vec3( 0.0, -1000.0, 0.0 );
    AddOption("translation","how much to move hanging strands by per timestep", Vec3( -0.05, 0.0, 0.0 ) );

    // Rod opts
    GetIntOpt("nv") = 10;
    
    // Additional parameters:
    AddOption("totalLength","enforced strand length", 20.0 );
    
    // Pre-setup to default values:
    GetScalarOpt( "stochasticPruningFraction" ) = 0.5;
    GetBoolOpt("useProxRodRodCollisions") = true;
    GetScalarOpt("collisionRadius") = 0.15; 
    GetScalarOpt( "externalCollisionsRadius" ) = 5.5;

    GetBoolOpt("useCTRodRodCollisions") = true;
    
    GetIntOpt("numberOfThreads") = 3;
    GetScalarOpt("strand_mu") = 0.3;
    AddOption("time_moving","Run along rods for the first t secs", 0.16 );
}

Aleka::~Aleka()
{}

std::vector<bool> frozen;
void Aleka::setupStrands()
{
    // discrete rod params
    const int nVertices = GetIntOpt("nv");
    const int nDOFs = 4 * nVertices - 1;
    
    // rod params
    const Scalar totalLength = GetScalarOpt("totalLength");
    const Scalar radiusA = GetScalarOpt("radius");
    const Scalar youngsModulus = GetScalarOpt("youngs-modulus");
    const Scalar shearModulus = GetScalarOpt("shear-modulus");
    const Scalar density = GetScalarOpt("density");
    const Scalar airDrag = GetScalarOpt("air-drag");
    const Scalar viscosity = GetScalarOpt("viscosity");
    const Scalar baseRotation = GetScalarOpt("base-rotation");
    
    const int fixed_rod_count = 5;
    AddOption("fixed_rod_count","", fixed_rod_count );
    const int moving_rod_count = 10;
    const Scalar fixed_rod_offset = 0.3;
    const Scalar moving_rod_offset = 0.3;
    
    const Scalar fixed_width = fixed_rod_count * ( 2 * radiusA + fixed_rod_offset );
    const Scalar fixed_inc = fixed_width / fixed_rod_count;
    const Scalar moving_width = moving_rod_count * ( 2 * radiusA + moving_rod_offset );
    const Scalar moving_inc = moving_width / moving_rod_count;
    const Scalar moving_x_offset = (fixed_width / 2) + 2;
    const Scalar moving_y_offset = 7.5;
    const Scalar moving_z_offset = 0;
    
    int layerTotal = fixed_rod_count + moving_rod_count;
    int layers = 1;

    frozen.resize( layerTotal );
    for( int layer = 0; layer < layers; ++layer ){

    int rod_id;
    // make fixed rods
    for ( rod_id = 0; rod_id < fixed_rod_count; ++rod_id )
    {
        // Prepare initial rod/strand position
        Vec3Array i_vertices;
        
        // Store arbitrary vertex coordinates
        for ( int i = 0; i < nVertices; ++i )
        {
            i_vertices.push_back( Vec3d( fixed_inc * rod_id - (fixed_width / 2),
                                         layer * (fixed_inc * 1) + (fixed_rod_count - rod_id) * fixed_inc,
                                         totalLength * i / ( nVertices - 1 ) - (totalLength/2) ) );
        }
        
        // Enforce the total length
        Scalar length = 0.0;
        for ( int i = 0; i < nVertices - 1; i++ )
            length += ( i_vertices[i + 1] - i_vertices[i] ).norm();
        for ( int i = 0; i < nVertices; i++ )
            i_vertices[i] *= totalLength / length;
        
        VecXd dofs( nDOFs );
        
        // Initial strand position
        for ( int i = 0; i < dofs.size(); i += 4 )
            dofs.segment<3>( i ) = i_vertices[ i / 4 ];
        
        DOFScriptingController* controller = new DOFScriptingController( );
        controller->freezeVertices( 0, true );
        controller->freezeVertices( nVertices - 1 );
        
        frozen[ layer * layerTotal + rod_id ] = true;

        ElasticStrandParameters* params = new ElasticStrandParameters( radiusA, youngsModulus, shearModulus, density, viscosity, airDrag, baseRotation );
        
        ElasticStrand* strand = new ElasticStrand( dofs, *params, controller );
        strand->setGlobalIndex( layer * layerTotal + rod_id );
        setRodCollisionParameters( *strand );
        m_strands.push_back( strand );
        
    }

    // make moving rods
    for( ; rod_id < moving_rod_count + fixed_rod_count; ++rod_id )
    {

        // Prepare initial rod/strand position
        Vec3Array i_vertices;
        
        // Store arbitrary vertex coordinates
        for ( int i = 0; i < nVertices; ++i )
        {
//            //Flat --> falling
//            i_vertices.push_back( Vec3d( moving_x_offset + totalLength * i / ( nVertices - 1 ) - (totalLength/2),
//                                         moving_y_offset,
//                                         moving_z_offset + moving_inc * ( rod_id - fixed_rod_count ) - (moving_width / 2) ) );
            //Hanging --> scripted
            i_vertices.push_back( Vec3d( moving_x_offset + layer * (moving_inc * 2),
                                         moving_y_offset - totalLength * i / ( nVertices - 1 ),
                                         moving_z_offset + moving_inc * ( rod_id - fixed_rod_count ) - (moving_width / 2) ) );
        }
        
        // Enforce the total length
        Scalar length = 0.0;
        for ( int i = 0; i < nVertices - 1; i++ )
            length += ( i_vertices[i + 1] - i_vertices[i] ).norm();
        for ( int i = 0; i < nVertices; i++ )
            i_vertices[i] *= totalLength / length;
        
        VecXd dofs( nDOFs );
        
        // Initial strand position
        for ( int i = 0; i < dofs.size(); i += 4 )
            dofs.segment<3>( i ) = i_vertices[ i / 4 ];
        
        DOFScriptingController* controller = new DOFScriptingController();
        controller->freezeVertices( 0, true );
        
        ElasticStrandParameters* params = new ElasticStrandParameters( radiusA, youngsModulus, shearModulus, density, viscosity, airDrag, baseRotation );
        
        frozen[ layer * layerTotal + rod_id ] = false;


        ElasticStrand* strand = new ElasticStrand( dofs, *params, controller );
        strand->setGlobalIndex( layer * layerTotal + rod_id );
        setRodCollisionParameters( *strand );
        m_strands.push_back( strand );

    }

    }
    
    std::cout << "num strands = " << m_strands.size() << std::endl;
    std::cout << "num dofs per strand = " << nDOFs << std::endl << std::endl;
    
    // global body forces
    GravitationForce::setGravity( GetVecOpt("gravity").cast<Scalar>() );
}

void Aleka::setupMeshes()
{}

int frame = 0;
bool Aleka::executeScript()
{
    if( getTime() > 8 ){
        std::exit(0);
    }

    Vec3 zero(  0. , 0., 0. );
    Mat3x id = Mat3x::Identity();

    const int fixed_rod_count = GetIntOpt("fixed_rod_count");
    int count = 0;
    for(auto rd_itr = m_strands.begin(); rd_itr != m_strands.end(); ++ rd_itr)
    {
        if( !frozen[count] ){

            // if ( getTime() < GetScalarOpt("time_moving") ){
            if ( getTime() < 0.25 ){
            // if ( (frame / 10 ) % 3 == 0 ){
                Vec3 translate = GetVecOpt("translation");
                SceneUtils::transformRodRootVtx( *rd_itr, id, zero, translate, 0 );
            }
            else{
                SceneUtils::transformRodRootVtx( *rd_itr, id, zero, zero, 0 );
            }
        }
        ++count;
    }

    ++frame;
    return true;
}

