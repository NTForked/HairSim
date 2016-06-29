#include "SingleContact.hh"

SingleContact::SingleContact() :
Scene("SingleContact", "Single Contact, simple orthogonal setup"),
m_radius(3.)
{
    AddOption( m_problemName, m_problemDesc, "" );
    
    // Global opts
    GetScalarOpt("dt") = 1e-3;
    GetVecOpt("gravity") = Vec3( 0.0, -100.0, 0.0 );
    // GetVecOpt("gravity") = Vec3( 0.0, 0.0, 0.0 );

    // Rod opts
    GetIntOpt("nv") = 6; //try 11, 981.0, uneven vertices (shows it better)
     
    // Additional parameters:
    AddOption("totalLength","enforced strand length", 10.0 );
    AddOption("fix_all_verts","script all vertices to be fixed", false );
    AddOption("x_offset","between 0 to totalLength/2", 0.0 );
    //  AddOption("thickness", "thickness", 0.5);
    //  AddOption("stiffness", "stiffness", 1000.0);
    //  AddOption("mass-damping", "mass damping for the rod", 0.0);

    // Pre-setup to default values:
    GetScalarOpt( "stochasticPruningFraction" ) = 0.5;
    GetBoolOpt("useProxRodRodCollisions") = true;
    GetScalarOpt("selfCollisionsRadius") = 0.16; //std::max(GetScalarOpt("radiusA"), GetScalarOpt("radiusB"));
    GetScalarOpt( "externalCollisionsRadius" ) = 0.75; /// thickness for CCD tests CULLING bounding boxes

    GetBoolOpt("useCTRodRodCollisions") = true;
    GetScalarOpt("percentCTRodRodCollisionsAccept") = 100.0;
    
    GetIntOpt("numberOfThreads") = 1;

    GetScalarOpt("strand_mu") = 0.3;
    setKeyBoardScriptingStrand( 1 );
}

SingleContact::~SingleContact()
{   
}

void SingleContact::setupStrands()
{
    std::cout << "Scene:: SingleContact" << std::endl;

    // discrete rod params
    const int nVertices = GetIntOpt("nv");
    const int nDOFs = 4 * nVertices - 1;
    
    // rod params
    const Scalar totalLength = GetScalarOpt("totalLength");
    const Scalar x_offset = GetScalarOpt("x_offset");
    const Scalar radiusA = GetScalarOpt("radiusA");
    const Scalar radiusB = GetScalarOpt("radiusB");
    const Scalar youngsModulus = GetScalarOpt("youngs-modulus");
    const Scalar shearModulus = GetScalarOpt("shear-modulus");
    const Scalar density = GetScalarOpt("density");
    const Scalar airDrag = GetScalarOpt("air-drag");
    const Scalar viscosity = GetScalarOpt("viscosity");
    const Scalar baseRotation = GetScalarOpt("base-rotation");
        
    // make rods
    for (int rod_id = 0; rod_id < 2; ++rod_id)
    {
        // Prepare initial rod/strand position
        Vec3Array i_vertices;
        
        // Store arbitrary vertex coordinates
        for ( int i = 0; i < nVertices; i++ )
        {
            if( rod_id == 0 ){
                i_vertices.push_back( Vec3d( x_offset, 0., totalLength * i / ( nVertices - 1 ) - (totalLength/2) ) );
            }
            else if( rod_id == 2 ){
                i_vertices.push_back( Vec3d( totalLength * i / ( nVertices - 1 ) - (totalLength/2), 1.5, 0.4 ) );
            }
            else{
                i_vertices.push_back( Vec3d( totalLength * i / ( nVertices - 1 ) - (totalLength/2), 0.5, 0. ) );
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
        
        Vec3Array scripted_vertices;
        scripted_vertices.push_back( i_vertices[0] );
        DOFScriptingController* controller = new DOFScriptingController( scripted_vertices );

        if( rod_id == 0 ){
            controller->freezeRootVertices<1>();            
            controller->freezeVertices( nVertices - 1 );
        }
        
        ElasticStrandParameters* params = new ElasticStrandParameters( 
                                                radiusA, 
                                                radiusB, 
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
        
        // extra stuff for render, etc...
        RodData* rd = new RodData( *strand, *controller );
        m_rodDatum.push_back( rd );
    }

    std::cout << "num strands = " << m_strands.size() <<'\n';
    std::cout << "num dofs per strand = " << nDOFs <<'\n';
    
    GravitationForce::setGravity( GetVecOpt("gravity").cast<Scalar>() );
}

void SingleContact::setupMeshes()
{
    // make collision mesh
    SimpleMeshController* mesh_controller = new SimpleMeshController( 0., m_dt );
    m_meshScripting_controllers.push_back( mesh_controller ); 
}

bool SingleContact::executeScript()
{

    // if( getTime() >= 0.5 ){
    //     GetVecOpt("gravity") = Vec3( 0.0, 0.0, 0.0 );
    // }


    // static int count = 1;
    // int vidx = 2;
    // Vec3 translate(0.0, 0.1, 0.0);
    // Vec3 zero(0.0, 0.0, 0.0);
    // for(auto rd_itr = m_rodDatum.begin(); rd_itr != m_rodDatum.end(); ++rd_itr)
    // {
    //     if( (*rd_itr)->getStrand().getGlobalIndex() == 1 ) continue;
    //     if ( count % 500 == 0 ){
    //         translateRodVertex( **rd_itr, vidx, translate );
    //         translateRodVertex( **rd_itr, vidx + 1, translate );
    //     }
    //     else{
    //         (*rd_itr)->getDofController().m_scriptedDegreesOfFreedom.erase(4*vidx + 0);            
    //         (*rd_itr)->getDofController().m_scriptedDegreesOfFreedom.erase(4*vidx + 1);            
    //         (*rd_itr)->getDofController().m_scriptedDegreesOfFreedom.erase(4*vidx + 2);  
    //         (*rd_itr)->getDofController().m_scriptedDegreesOfFreedom.erase(4*3 + 0);            
    //         (*rd_itr)->getDofController().m_scriptedDegreesOfFreedom.erase(4*3 + 1);            
    //         (*rd_itr)->getDofController().m_scriptedDegreesOfFreedom.erase(4*3 + 2);                        
    //     }
    //         // unscriptRodVertex( **rd_itr, vidx );
    // }

    // ++count;
    return true;
}

