#include "MultipleContact.h"

MultipleContact::MultipleContact() :
Scene("MultipleContact", "MultipleContact, one strand edge falls between fixed strands with opposing normals"),
m_radius(3.)
{
    AddOption( m_problemName, m_problemDesc, "" );

    // Global opts
    GetScalarOpt("dt") = 1e-3;
    GetVecOpt("gravity") = Vec3( 0.0, -1000.0, 0.0 );

    // Rod opts
    GetIntOpt("nv") = 10; //try 11 for uneven vertices (Corner Case)
    
    // Additional parameters:
    AddOption("totalLength","enforced strand length", 20.0 );
    AddOption("fix_all_verts","script all vertices to be fixed", false );
    AddOption("x_offset","between 0 to totalLength/2", 0.0 );
    //  AddOption("thickness", "thickness", 0.5);
    //  AddOption("stiffness", "stiffness", 1000.0);
    //  AddOption("mass-damping", "mass damping for the rod", 0.0);

    // Pre-setup to default values:
    GetScalarOpt( "stochasticPruningFraction" ) = 0.5;
    GetBoolOpt("useProxRodRodCollisions") = true;
    GetScalarOpt("collisionRadius") = 0.16;
    GetScalarOpt( "externalCollisionsRadius" ) = 5.5;

    GetBoolOpt("useCTRodRodCollisions") = true;
    
    GetIntOpt("numberOfThreads") = 1; // 5
    GetScalarOpt("strand_mu") = 0.0;
}

MultipleContact::~MultipleContact()
{   
}

void MultipleContact::setupStrands()
{
    std::cout << "Scene:: MultipleContact" << std::endl;

    // discrete rod params
    const int nVertices = GetIntOpt("nv");
    const int nDOFs = 4 * nVertices - 1;
    
    // rod params
    const Scalar totalLength = GetScalarOpt("totalLength");
    const Scalar x_offset = GetScalarOpt("x_offset");
    const Scalar radiusA = GetScalarOpt("radius");
    const Scalar youngsModulus = GetScalarOpt("youngs-modulus");
    const Scalar shearModulus = GetScalarOpt("shear-modulus");
    const Scalar density = GetScalarOpt("density");
    const Scalar airDrag = GetScalarOpt("air-drag");
    const Scalar viscosity = GetScalarOpt("viscosity");
    const Scalar baseRotation = GetScalarOpt("base-rotation");
    
    bool flat = false;
    
    // make rods
    for (int rod_id = 0; rod_id < 4; ++rod_id)
    {
        // Prepare initial rod/strand position
        Vec3Array i_vertices;
        
        // Store arbitrary vertex coordinates
        for ( int i = 0; i < nVertices; i++ )
        {
            
            if( rod_id == 0 ){
                i_vertices.push_back( Vec3d( totalLength * i / ( nVertices - 1 ) - (totalLength/2), 5, 0. ) );
            }
            else if( flat ) {
                if ( rod_id == 1 ){
                    i_vertices.push_back( Vec3d( -0.6, 0.0, totalLength * i / ( nVertices - 1 ) - (totalLength/2) ) );
                }
                else if ( rod_id == 2 ){
                    i_vertices.push_back( Vec3d( 0.0, 0.0, ( totalLength * i / ( nVertices - 1 ) - (totalLength/2))));
                }
                else if ( rod_id == 3 ){
                    i_vertices.push_back( Vec3d( 0.6, 0.0, totalLength * i / ( nVertices - 1 ) - (totalLength/2) ) );
                }
            }
            else{
                if ( rod_id == 1 ){
                    i_vertices.push_back( Vec3d( -0.6, totalLength * i / ( nVertices - 1 ) - (totalLength/2), totalLength * i / ( nVertices - 1 ) - (totalLength/2) ) );
                }
                else if ( rod_id == 2 ){
                    i_vertices.push_back( Vec3d( 0., totalLength -  totalLength * i / ( nVertices - 1 ) - (totalLength/2), ( totalLength * i / ( nVertices - 1 ) - (totalLength/2))));
                }
                else if ( rod_id == 3 ){
                    i_vertices.push_back( Vec3d( 0.6, totalLength * i / ( nVertices - 1 ) - (totalLength/2), totalLength * i / ( nVertices - 1 ) - (totalLength/2) ) );
                }
            }
        }
        
        // Enforce the total length
        Scalar length = 0.0;
        for ( int i = 0; i < nVertices - 1; i++ )
            length += ( i_vertices[i + 1] - i_vertices[i] ).norm();
        std::cout << "actual_length: " << length <<  std::endl;
        for ( int i = 0; i < nVertices; i++ )
            i_vertices[i] *= totalLength / length;
        
        VecXd dofs( nDOFs );
        
        // Initial strand position
        for ( int i = 0; i < dofs.size(); i += 4 )
            dofs.segment<3>( i ) = i_vertices[i / 4];
        
        DOFScriptingController* controller = new DOFScriptingController( );
        
        if( rod_id != 0 ){
            controller->freezeVertices(0, true);

            if( GetBoolOpt("fix_all_verts") ){
                for( int b = 0; b < nVertices; ++b ){
                    controller->freezeVertices( b );
                }
            }
            
            controller->freezeVertices( nVertices - 1, true );
        }
        
        ElasticStrandParameters* params = new ElasticStrandParameters( 
                                                radiusA, 
                                                youngsModulus, 
                                                shearModulus, 
                                                density, 
                                                viscosity, 
                                                airDrag, 
                                                baseRotation );
        
        ElasticStrand* strand = new ElasticStrand( dofs, *params, controller );
        strand->setGlobalIndex( rod_id );
        setRodCollisionParameters( *strand );
        m_strands.push_back( strand );

    }

    std::cout << "num strands = " << m_strands.size() <<'\n';
    std::cout << "num dofs per strand = " << nDOFs <<'\n';
    
    GravitationForce::setGravity( GetVecOpt("gravity").cast<Scalar>() );
}

void MultipleContact::setupMeshes()
{

}

bool MultipleContact::executeScript()
{
    return true;
}

