#include "CousinIt.hh"

CousinIt::CousinIt() :
ProblemStepper("Cousin It", "spinning sphere")
{
    GetIntOpt("nv") = 25;
    AddOption("curl_radius", "", 0.2 );
    AddOption("curl_density", "", 2.2 );
    GetScalarOpt("curl_density_perturbation") = 0.0001;
    GetScalarOpt("curl_radius_perturbation") = 0.0001;
    GetScalarOpt("density") = 1.32;
    GetScalarOpt("youngs-modulus") = 3.9e+09;
    GetScalarOpt("viscosity") = 5e+07;
    GetScalarOpt("stretchingThreshold") = 1.;
    GetScalarOpt("radiusA") = 0.0037;
    GetScalarOpt("radiusB") = 0.0037;
    GetIntOpt("maxNewtonIterations") = 30;
    AddOption("geometryscale","uniform rescaling of sphere and hair root locations", 2.25);
    // GetBoolOpt("failure_testing_on") = true;


    AddOption("root_length", "", 0.1);
    AddOption("hair_length", "", 17.);
  
    AddOption("scripting_on","", true );
    AddOption("end_time","", 10.);
    
    AddOption("hair_region", "(0, 1] specifies percent of sphere to cover, 0 specifies load hair", 0.5 );
    AddOption("hair_load_file","if hair_region  = 0", "assets/problemfiles/cousinit_profile_startRods.ply");
    AddOption("num_hairs", "number of hair seeds to grow in hair_region (i.e., if hair_region != 0) ", 1000);
    AddOption("rotatescale","uniform rescaling of rotational speeds", 1.0);
    AddOption("sphere_mesh_filename","","assets/TriangulatedSphere.obj"); 
    
    AddOption("time_offset","start scripting at time t + offset", 0.);

    // GetBoolOpt("render") = false;
    GetIntOpt("numberOfThreads") = 3; // 32?

    // Pre-setup to default values:
    GetScalarOpt( "stochasticPruningFraction" ) = 0.5;
    GetBoolOpt("useProxRodRodCollisions") = true;
    GetBoolOpt("useCTRodRodCollisions") = false;
    GetScalarOpt("selfCollisionsRadius") = 0.0025;
    GetScalarOpt("percentCTRodRodCollisionsAccept") = 100.0;

    GetStringOpt("checkpointDir") = "/Users/henrique/Desktop/LearnHair/build/Apps/StrandSimulator";

    GetScalarOpt("dt") = 0.05; // 1e-2;
    GetScalarOpt("strand_mu") = 0.3;
    GetScalarOpt("mesh_mu") = 0.1;
}

CousinIt::~CousinIt()
{
    
}

void CousinIt::generateNormalSamples( Scalar hair_region, int num_hairs, std::vector<Vec3x>& normals )
{
    typedef boost::minstd_rand rng_type;
    typedef boost::uniform_real<double> dist_type;
    typedef boost::variate_generator<rng_type&, dist_type> generator_type;
    

    // hypercube rejection sampling
    dist_type dist1( 1 - 2 * hair_region, 1 );
    rng_type generator1;
    generator_type gen1(generator1, dist1);
    dist_type dist2( -1, 0.99); //1 );
    rng_type generator2;
    generator_type gen2(generator2, dist2);
    
    normals.clear();
    for( int j = 0; j < num_hairs; ++j )
    {
        Scalar y = gen2();//gen1();
        Scalar x = gen2();
        Scalar z = gen2();
        Scalar a = sqrt(x*x + y*y + z*z);
        while (a > 1 || y/a < 1 - 2 * hair_region)
        {
            y = gen2();//gen1();
            x = gen2();
            z = gen2();
            a = sqrt(x*x + y*y + z*z);
        }
        
        Vec3d newnormal(x/a,y/a,z/a);
        normals.push_back(newnormal);
    }
}


void CousinIt::loadStrands()
{
    std::string hair_config_filename = GetStringOpt("hair_load_file");
    std::stringstream name;
    name << hair_config_filename;
    std::ifstream infile(name.str().c_str());
    
    if( !infile.is_open() )
    {
        std::cerr << "\033[31;1mCousinIt:\033[m can't load " << hair_config_filename << std::endl;
        exit(1);
    }
    
    //jump over ply header
    std::string header_string;
    for(int i=0; i< 11; i++) getline( infile, header_string);
    
    int rod_id = 0;
    int current_segment = 0;
    int segment = current_segment;
    double x, y, z;
    while(!infile.eof())
    {
        if ( !( rod_id < GetIntOpt("num_hairs") ) )
            break;
        
        std::vector<Vec3x> vertices;
        if (segment != current_segment)
        {
            Vec3x pos(x,y,z);
            vertices.push_back(pos);
            current_segment = segment;
        }
        while(true)
        {
            std::string vert_string;
            getline( infile, vert_string);
            std::istringstream vertstream(vert_string);
            vertstream >> x >> y >> z >> segment;
            if(segment == current_segment && !infile.eof())
            {
                Vec3x pos(x,y,z);
                vertices.push_back(pos);
            }
            else{
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
//                // initial strand position
//                int num_DoFs = 4 * vertices.size() - 1;
//                VecXd dofs( num_DoFs );
//                for ( int i = 0; i < dofs.size(); i += 4 )
//                {
//                    dofs.segment<3>( i ) = vertices[i / 4];
//                    if( i+3 < dofs.size() ) dofs( i + 3 ) = 0.;
//                }
                
                strandsim::Vec3xArray scripted_vertices;
                scripted_vertices.push_back( vertices[0] );
                scripted_vertices.push_back( vertices[1] );
                DOFScriptingController* controller = new DOFScriptingController( scripted_vertices );
                controller->freezeRootVertices<2>();
                
                ElasticStrandParameters* params = new ElasticStrandParameters( radiusA, radiusB, youngsModulus, shearModulus, density, viscosity, airDrag, baseRotation );
                
                ElasticStrand* strand = new ElasticStrand( dofs, *params, controller );
                strand->setGlobalIndex( rod_id );
                setRodCollisionParameters( *strand );
                m_strands.push_back( strand );
                
                // extra stuff for render, etc...
                RodData* rd = new RodData( *strand, *controller );
                m_rodDatum.push_back( rd );
                ++rod_id;
                
                //std::cout<< "loaded rod with dofs = "<< dofs << std::endl;
                
                break;
            }
        }
    }
    
    std::cout << "# \033[35;1mCousinIt message:\033[m loaded " << m_rodDatum.size() << " rods." << std::endl;
}



void CousinIt::setupStrands()
{
    std::cout << "NEED to use Bogus, with and without friction, multi-threaded, loop on, no redundancies, real abscissas, render should be off, dump coords appropriately, vertices displayed off" << std::endl;
    if( GetScalarOpt("hair_region") == 0. )
        loadStrands();
    else {
        // rod options
        Scalar radiusA = GetScalarOpt("radiusA");
        Scalar radiusB = GetScalarOpt("radiusB");
        Scalar youngsModulus = GetScalarOpt("youngs-modulus");
        Scalar shearModulus = GetScalarOpt("shear-modulus");
        Scalar density = GetScalarOpt("density");
        Scalar airDrag = GetScalarOpt("air-drag");
        Scalar viscosity = GetScalarOpt("viscosity");
        Scalar baseRotation = GetScalarOpt("base-rotation");
        //  Scalar thickness = GetScalarOpt("thickness");
        //  Scalar stiffness = GetScalarOpt("stiffness");
        
        // sample rod positions
        std::vector<Vec3x> initialnormals;
        generateNormalSamples( GetScalarOpt("hair_region"), GetIntOpt("num_hairs"), initialnormals );
        
        // sphere geometry
        double geo_rescale = GetScalarOpt("geometryscale");
        double sphereradius = geo_rescale * 4.0; //cm
        
        // more rod params
        int num_vertices = GetIntOpt("nv");
        int num_DoFs = 4 * num_vertices - 1;
        double L = GetScalarOpt("hair_length");
        double dL = L / ( (double) num_vertices - 1 );
        
        int rod_id = 0;
        for( int i = 0; i < (int) initialnormals.size(); ++i )
        {
            std::vector<Vec3x> vertices;
            
            Scalar curl_radius = GetScalarOpt("curl_radius");
            Scalar curl_density = GetScalarOpt("curl_density");
            Scalar root_length = GetScalarOpt("root_length");
            
            if( curl_radius == 0. || curl_density == 0. )
                generateStraightHair( initialnormals[i], sphereradius * initialnormals[i], dL, num_vertices, vertices, root_length );
            else
                generateCurlyHair( initialnormals[i], sphereradius * initialnormals[i], dL, num_vertices, vertices, curl_radius, curl_density, root_length );
            
            // initial strand position
            VecXd dofs( num_DoFs );
            for ( int i = 0; i < dofs.size(); i += 4 )
                dofs.segment<3>( i ) = vertices[i / 4];
            
            strandsim::Vec3xArray scripted_vertices;
            scripted_vertices.push_back( vertices[0] );
            scripted_vertices.push_back( vertices[1] );
            DOFScriptingController* controller = new DOFScriptingController( scripted_vertices );
            controller->freezeRootVertices<2>();
            
            ElasticStrandParameters* params = new ElasticStrandParameters( radiusA, radiusB, youngsModulus, shearModulus, density, viscosity, airDrag, baseRotation );
            
            ElasticStrand* strand = new ElasticStrand( dofs, *params, controller, GetScalarOpt("selfCollisionsRadius") );
            strand->setGlobalIndex( rod_id );
            setRodCollisionParameters( *strand );
            m_strands.push_back( strand );
            
            // extra stuff for render, etc...
            RodData* rd = new RodData( *strand, *controller );
            m_rodDatum.push_back( rd );
            ++rod_id;
        }
        std::cout << "# num strands: \n" << m_strands.size() <<'\n';
        std::cout << "# num dofs per strand: \n" << num_DoFs <<'\n';
    }
    
}

void CousinIt::setupMeshes()
{
    // make collision mesh
    SimpleMeshController* mesh_controller = new SimpleMeshController( 0., m_dt );
    m_meshScripting_controllers.push_back( mesh_controller );
    std::string file_name = GetStringOpt( "sphere_mesh_filename" ) ;
    mesh_controller->loadMesh(file_name);
    // set mu for mesh
    mesh_controller->setDefaultFrictionCoefficient( GetScalarOpt( "mesh_mu" ) ) ;
    
    // rescale sphere -- assumes mesh is centered at origin
    TriangularMesh* currentMesh = m_meshScripting_controllers[0]->getCurrentMesh();
    double meshradius = (currentMesh->getVertex(0)).norm();
    //std::cout<< "radius before = " << meshradius << '\n';
    double geo_rescale = GetScalarOpt("geometryscale");
    double sphereradius = geo_rescale * 4.0;
    double meshscale = sphereradius/meshradius;
    //std::cout<< "mesh scale = " << meshscale << '\n';
    for( unsigned i = 0; i < currentMesh->nv(); ++i )
    {
        currentMesh->setVertex(i, currentMesh->getVertex(i) * meshscale);
    }
    meshradius = (currentMesh->getVertex(0)).norm();
    std::cout<< "# sphere radius: \n" << meshradius << '\n';
    
    strandsim::TriangleMeshRenderer* mesh_renderer = new strandsim::TriangleMeshRenderer( *(mesh_controller->getCurrentMesh()));
    m_mesh_renderers.push_back( mesh_renderer );
}

bool CousinIt::executeScript()
{
    Vec3x zero(0.,0.,0.);
    
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
  
    // std::cout << " scripting_on is " << GetBoolOpt("scripting_on") << "\n";
  
    if ( GetBoolOpt("scripting_on") )
    {
      if( getTime() <= 2. - GetScalarOpt("time_offset") )
      {
        freezeTriangleObject( *currentMesh );
        for(auto rd_itr = m_rodDatum.begin(); rd_itr != m_rodDatum.end(); ++rd_itr)
        {
          freezeRodRoot( **rd_itr );
        }
      }
      else if ( getTime() <= 3 - GetScalarOpt("time_offset") )
      {
        // transform mesh
        transformTriangleObject( *currentMesh, rotz, zero, zero );
        // transform rod BC
        for(auto rd_itr = m_rodDatum.begin(); rd_itr != m_rodDatum.end(); ++rd_itr)
        {
          transformRodRoot( **rd_itr, rotz, zero, zero );
        }
      }
      else if(getTime() <= 4. - GetScalarOpt("time_offset") )
      {
        freezeTriangleObject( *currentMesh );
        for(auto rd_itr = m_rodDatum.begin(); rd_itr != m_rodDatum.end(); ++rd_itr)
        {
          freezeRodRoot( **rd_itr );
        }
      }
      else if(getTime() <= 5 - GetScalarOpt("time_offset") )
      {
        // back the other way
        // transform mesh
        transformTriangleObject( *currentMesh, rotzT, zero, zero );
        // transform rod BC
        for(auto rd_itr = m_rodDatum.begin(); rd_itr != m_rodDatum.end(); ++rd_itr)
        {
          transformRodRoot( **rd_itr, rotzT, zero, zero );
        }
      }
      else if(getTime() <= 6. - GetScalarOpt("time_offset") )
      {
        freezeTriangleObject( *currentMesh );
        for(auto rd_itr = m_rodDatum.begin(); rd_itr != m_rodDatum.end(); ++rd_itr)
        {
          freezeRodRoot( **rd_itr );
        }
      }
      else if(getTime() <= 7 - GetScalarOpt("time_offset") )
      {
        // transform mesh
        transformTriangleObject( *currentMesh, rotx, zero, zero );
        // transform rod BC
        for(auto rd_itr = m_rodDatum.begin(); rd_itr != m_rodDatum.end(); ++rd_itr)
        {
          transformRodRoot( **rd_itr, rotx, zero, zero );
        }
      }
      else if(getTime() <= 8. - GetScalarOpt("time_offset") )
      {
        freezeTriangleObject( *currentMesh );
        for(auto rd_itr = m_rodDatum.begin(); rd_itr != m_rodDatum.end(); ++rd_itr)
        {
          freezeRodRoot( **rd_itr );
        }
      }
      else if(getTime() <= 9 - GetScalarOpt("time_offset") )
      {
        // back the other way
        // transform mesh
        transformTriangleObject( *currentMesh, rotxT, zero, zero );
        // transform rod BC
        for(auto rd_itr = m_rodDatum.begin(); rd_itr != m_rodDatum.end(); ++rd_itr)
        {
          transformRodRoot( **rd_itr, rotxT, zero, zero );
        }
      }
      else if(getTime() <= 10. - GetScalarOpt("time_offset") )
      {
        freezeTriangleObject( *currentMesh );
        for(auto rd_itr = m_rodDatum.begin(); rd_itr != m_rodDatum.end(); ++rd_itr)
        {
          freezeRodRoot( **rd_itr );
        }
      }
      else if(getTime() <= 11 - GetScalarOpt("time_offset") )
      {
        // transform mesh
        transformTriangleObject( *currentMesh, roty, zero, zero );
        // transform rod BC
        for(auto rd_itr = m_rodDatum.begin(); rd_itr != m_rodDatum.end(); ++rd_itr)
        {
          transformRodRoot( **rd_itr, roty, zero, zero );
        }
      }
      else if(getTime() <= 12. - GetScalarOpt("time_offset")  )
      {
        freezeTriangleObject( *currentMesh );
        for(auto rd_itr = m_rodDatum.begin(); rd_itr != m_rodDatum.end(); ++rd_itr)
        {
          freezeRodRoot( **rd_itr );
        }
      }
      else if(getTime() <= 13 - GetScalarOpt("time_offset") )
      {
        // back the other way
        // transform mesh
        transformTriangleObject( *currentMesh, rotyT, zero, zero );
        // transform rod BC
        for(auto rd_itr = m_rodDatum.begin(); rd_itr != m_rodDatum.end(); ++rd_itr)
        {
          transformRodRoot( **rd_itr, rotyT, zero, zero );
        }
      }
      else
      {
        freezeTriangleObject( *currentMesh );
        for(auto rd_itr = m_rodDatum.begin(); rd_itr != m_rodDatum.end(); ++rd_itr)
        {
          freezeRodRoot( **rd_itr );
        }
      }
    }
    else{
      freezeTriangleObject( *currentMesh );
      for(auto rd_itr = m_rodDatum.begin(); rd_itr != m_rodDatum.end(); ++rd_itr)
      {
        freezeRodRoot( **rd_itr );
      }
    }

    if ( getTime() >= GetScalarOpt("end_time") )
    {
        std::cout << "# Simulation complete. Exiting." << std::endl;
        exit(0);
    }
  
    return true;
  
}

