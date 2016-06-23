#include "PlayBack.hh"

PlayBack::PlayBack() :
ProblemStepper("Play Back", "play back output files"),
m_current_frame(0)
{

AddOption("num_frames", "number of frames to play back", 10001); 
AddOption("hair_thick", "rendering thickness", 0.02); 
AddOption("dt_time", "dt used in original sim", 0.001); 

AddOption("padded_width", "number of 0s in frame files", 8); 

AddOption("realtime" , "skip set number of frames" , true );

// AddOption( "directory" , "file directory with ordered rod (single ply) and mesh (numbered obj) files" , "/Users/henrique/Desktop/LearnHair/build/Apps/StrandSimulator/2016_06_03_14_42_00_simulation_capture" );
AddOption( "directory" , "file directory with ordered rod (single ply) and mesh (numbered obj) files" , "/Users/henrique/Desktop/tmp/tmp" );
AddOption( "mesh_saved_elsewhere" , "if mesh is saved elsewhere" , false );
AddOption( "mesh_directory" , "when mesh_saved_elsewhere = true, file directory with ordered mesh (numbered obj) files" , "assets/TriangulatedSphere.obj" );
AddOption( "mesh_start_name" , "when mesh_saved_elsewhere = true, name of (numbered obj) files" , "TriangulatedSphere.obj" );
AddOption("start_frame", "number of frames to play back", 0);
AddOption("number_of_meshes","number of meshes to load", 0);
AddOption("number_of_rods","number of rods to load, -1 defaults to all", -1);
}

PlayBack::~PlayBack(){}

void PlayBack::loadRods(int frame)
{
    std::string inputdirectory = GetStringOpt("directory");
    
    int file_width = GetIntOpt("padded_width");
    std::stringstream name;
    name << inputdirectory << "/rods_" << std::setfill('0') << std::setw(file_width) << frame << ".ply";
    std::ifstream infile(name.str().c_str());
    
    if( !infile.is_open() )
    {
        std::cerr << "\033[31;1mPlayback:\033[m can't load " << name << std::endl;
        exit(1);
    }
    
    //jump over ply header
    std::string header_string;
    for(int i=0; i< 11; i++) getline( infile, header_string);
    
    if (frame == GetIntOpt("start_frame") )
    {
        int rod_id = 0;
        int current_segment = 0;
        int segment = current_segment;
        double x, y, z;
        while(!infile.eof())
        {
            if ( !( rod_id < GetIntOpt("number_of_rods") ) && GetIntOpt("number_of_rods") != -1 )
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
                    
                    strandsim::Vec3xArray scripted_vertices;
                    scripted_vertices.push_back( vertices[0] );
                    scripted_vertices.push_back( vertices[1] );
                    DOFScriptingController* controller = new DOFScriptingController( scripted_vertices );
                    controller->freezeRootVertices<2>();
                    
                    ElasticStrandParameters* params = new ElasticStrandParameters( radiusA, radiusB, youngsModulus, shearModulus, density, viscosity, airDrag, baseRotation );
                    
                    ElasticStrand* strand = new ElasticStrand( dofs, *params, controller, GetScalarOpt("hair_thick") );
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
    }
    else {
        int rod_id = 0;
        int current_segment = 0;
        int segment = current_segment;
        double x, y, z;
        while(!infile.eof())
        {
            if ( !( rod_id < GetIntOpt("number_of_rods") ) && GetIntOpt("number_of_rods") != -1 )
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
                    
                    for (int i = 0; i< (int) vertices.size(); ++i)
                    {
                        m_rodDatum[rod_id]->getStrand().setVertex( i, vertices[ i ] );
                    }
                    ++rod_id;
                    break;
                }
            }
        }
    }
    
    // std::cout << "\033[35;1mPlayBack message:\033[m loaded " << m_rodDatum.size() << " rods. from: " << name.str() << std::endl;
    std::cout << "\033[35;1mPlayBack message:\033[m loaded " << m_rodDatum.size() << " rods " << std::endl;
}

void PlayBack::loadMeshes(int frame)
{
    for(int mesh_num = 0; mesh_num < GetIntOpt("number_of_meshes"); ++mesh_num )
    {
        std::stringstream name;
        if( ! GetBoolOpt( "mesh_saved_elsewhere" ) )
        {
            std::string inputdirectory = GetStringOpt("directory");
            int file_width = 20;
            name << std::setfill('0');
            name << inputdirectory << "/mesh"<< mesh_num << "_" << std::setw(file_width) << frame << ".obj";
            std::cout<< "\033[31;1mPlayback:\033[m loading " << name.str() << std::endl;
        }
        else {
            std::string inputdirectory = GetStringOpt("mesh_directory");
            name << inputdirectory << "/"<< GetStringOpt("mesh_start_name") << "_" << ( frame + 2 )<< ".obj";
            std::cout<< "\033[31;1mPlayback:\033[m loading " << name.str() << std::endl;
        }
      
        std::string mesh_name = name.str();
        
        std::ifstream infile(name.str().c_str());
        
        if( !infile.is_open() )
        {
            std::cerr << "\033[31;1mPlayback:\033[m can't load " << name << std::endl;
            exit(1);
        }
        
        ObjParser objparser;
        if(frame == GetIntOpt("start_frame"))
        {
            SimpleMeshController* mesh = new SimpleMeshController( 0., m_dt );
            m_meshScripting_controllers.push_back( mesh );
            mesh->loadMesh( mesh_name );
            strandsim::TriangleMeshRenderer* mesh_renderer = new strandsim::TriangleMeshRenderer( *(mesh->getCurrentMesh()) );
            m_mesh_renderers.push_back( mesh_renderer );
        }
        else {
            //get next mesh
            SimpleMeshController* mesh = new SimpleMeshController( 0., m_dt );
            mesh->loadMesh( mesh_name );
            
            // set current mesh
            TriangularMesh* currentMesh = m_meshScripting_controllers[mesh_num]->getCurrentMesh();
            for( unsigned i = 0; i < currentMesh->nv(); ++i )
            {
                currentMesh->setVertex( i, mesh->getCurrentMesh()->getVertex(i) );
            }
            delete mesh;
            mesh = NULL;
        }
    }
    
    if( GetIntOpt("number_of_meshes") > 0 ){
        std::cout << "\033[35;1mPlayback message:\033[m loaded " <<  m_meshScripting_controllers.size() << " mesh(es)." << std::endl;
    }
}

void PlayBack::setupStrands()
{    
    // ensure screen grab outputs on every frame
    m_dt = 1. / GetScalarOpt("fps");

// Hairy ball with dt = 1e-2
//     m_dt = 1. / 50.0;
// int precision = 2;
// unsigned int uiTemp = (m_dt*pow(10, precision)) * 1; //Note I'm using unsigned int so that I can increase the precision of the truncate
// m_dt = (((double)uiTemp)/pow(10,precision) * 1);
    
    m_dt = GetScalarOpt("dt_time");

    m_current_frame = GetIntOpt("start_frame");
    std::cout << "\033[35;1mPlayback message:\033[m starting with frame " << GetIntOpt("start_frame") << std::endl;
    
    loadRods(m_current_frame);
    loadMeshes(m_current_frame);
    
    ++m_current_frame;
}

void PlayBack::setupMeshes()
{
    return;
    // make collision mesh
    SimpleMeshController* mesh_controller = new SimpleMeshController( 0., m_dt );
    m_meshScripting_controllers.push_back( mesh_controller );
    std::string file_name = GetStringOpt( "mesh_directory" ) ;
    mesh_controller->loadMesh(file_name);
    
    // rescale sphere -- assumes mesh is centered at origin
    TriangularMesh* currentMesh = m_meshScripting_controllers[0]->getCurrentMesh();
    double meshradius = (currentMesh->getVertex(0)).norm();
    //std::cout<< "radius before = " << meshradius << '\n';
    double geo_rescale = 2.25;
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

bool PlayBack::executeScript()
{

    if( GetBoolOpt("realtime") ){
        m_current_frame += 40; // 30 FPS playback
    }

    if( m_current_frame - GetIntOpt("start_frame") < GetIntOpt("num_frames") )
    {
        std::cout << "m_current_frame: " << m_current_frame << std::endl;
        loadRods(m_current_frame);
        loadMeshes(m_current_frame);
        ++m_current_frame;
    }
    else {
        {
            std::cout << "Playback complete. Exiting." << std::endl;
            exit(0);
        }
    }
    return true;
}




