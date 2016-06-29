#include "Scene.h"

#include <sys/stat.h>
#include <time.h>
#include <ios>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>

#include "../../StrandSim/Collision/ElementProxy.hh"
#include "../../StrandSim/Collision/CollisionDetector.hh"

Scene::ProblemStepper( const std::string& name, const std::string& desc )
: m_problemName( name )
, m_problemDesc( desc )
, m_t( 0. ) // start time
, m_isSimulated( true )
, m_strandsManager( NULL )
{
    addOptions();
    m_renderer = new StrandRenderer();
}

Scene::~ProblemStepper()
{
  if( m_strandsManager != NULL )
  {
    delete m_strandsManager;
    m_strandsManager = NULL;
  }

  for( std::vector< ElasticStrand* >::size_type i = 0; i < m_strands.size(); ++i )
  {
    if( m_strands[i] != NULL )
    {
      delete m_strands[i];
      m_strands[i] = NULL;
    }
  }
  
  for( std::vector< TriMesh* >::size_type i = 0; i < m_meshes.size(); ++i )
  {
    if( m_meshes[i] != NULL )
    {
      delete m_meshes[i];
      m_meshes[i] = NULL;
    }
  }
}

void Scene::setup()
{
    // stepping params :
    m_dt = GetScalarOpt( "dt" ); // ~frame
    std::cout << "# step_size: "<< m_dt << std::endl;
    
    setupStrands();
    setupMeshes();
    setSimulationParameters();
    
    if( m_isSimulated ){
        // Enforce desired or maximum number of threads
        {
            const int numThreads = m_simulation_params.m_numberOfThreads > 0 ? m_simulation_params.m_numberOfThreads : sysconf( _SC_NPROCESSORS_ONLN );
            omp_set_num_threads( numThreads );
        }
        m_strandsManager = new StrandImplicitManager( m_strands, m_simulation_params, m_meshes );
    }
}

bool Scene::step()
{
    std::cout << "\n# Stepping @ time: " << m_t << "\n";

    executeScript();
    if( m_isSimulated ){
        m_strandsManager->step( m_dt );
    }

    m_t += m_dt;    
    return true;
}

void Scene::clearContacts()
{
    m_renderer->verts.clear();
    m_renderer->displaceVec.setZero();
    m_renderer->arrows.clear();
    m_renderer->startProxies.clear();
    m_renderer->endProxies.clear(); 
    m_renderer->vertLabels.clear(); 
    m_renderer->twistLabels.clear(); 
}

void Scene::addOptions()
{ // Default options
    AddOption("checkpointDir", "directory where checkpointing file(s) stored", "");

    // world opts
    AddOption("dt", "time-step size", 1e-2);
    AddOption("gravity", "gravity", Vec3d(0, -981, 0));
    
    // rod opts
    AddOption("externalCollisionsRadius", "rod-mesh and rod-rod CT collision radius", 0.01 );
    AddOption("collisionRadius","rod-rod proximity collision radius", 0.01 );

    AddOption("radius", "radius of physical cross-section", 0.0025);
    AddOption("youngs-modulus", "Young's modulus for the rod", 1e10 );
    AddOption("shear-modulus", "Shear modulus for the rod", 3.4e9 );
    AddOption("density", "volumetric density of the rod", 1.3 );
    AddOption("air-drag", "", 1.e-4);
    AddOption("viscosity", "Dynamic viscosity for the rod", 1e3 );
    AddOption("base-rotation", "base-rotation", 0.);

    AddOption("nv", "number of vertices in the rod", 10 );
    
    AddOption( "curl_radius_perturbation" , "" , 5e-2 );
    AddOption( "curl_density_perturbation" , "" , 5e-2 );
    
    AddOption( "strand_mu" , "" , 0.3 );
    AddOption( "mesh_mu" , "" , 0.3 );

    AddOption( "rootImmunityLength" , "how much of the base of the strand to ignore" , 0.3 );
        
    // for implicit solve :
    AddOption("useGraphSplitOption", "whether to use graph split for large problems", false );
    AddOption("useLengthProjection", "whether to use length projections failsafe", false );
        
    AddOption("stretchingThreshold","exponent for stretch threshold", 2.); //2;
    AddOption("costretch_residual_threshold","exponent for costretch residual threshold", 0.);

    
    AddOption( "useProxRodRodCollisions" , "", true);
    AddOption( "useCTRodRodCollisions" , "", false );
    AddOption( "percentCTRodRodCollisionsAccept", "", 100. );
    AddOption( "useNonLinearAsFailsafe","", false );
    AddOption( "alwaysUseNonLinear","", true );

    AddOption("maxNewtonIterations","",10);
    
    // sys options
    AddOption("numberOfThreads","",4);
    AddOption("simulationManager_limitedMemory","", false);
    AddOption("contactProblem_limitedMemory","", false);
    
    //    // for YacFS solve :
    AddOption("gaussSeidelTolerance","", 1e-6 );
    AddOption("stochasticPruningFraction","factor for how many collisions to keep?", 0.8);
}

void Scene::setSimulationParameters()
{
    m_simulation_params.m_numberOfThreads = GetIntOpt("numberOfThreads");
    m_simulation_params.m_stretching_threshold = GetScalarOpt( "stretchingThreshold" );
    m_simulation_params.m_costretch_residual_threshold = GetScalarOpt( "costretch_residual_threshold" );
    m_simulation_params.m_useGraphSplitOption = GetBoolOpt( "useGraphSplitOption" );
    m_simulation_params.m_useLengthProjection = GetBoolOpt( "useLengthProjection" );
        
    m_simulation_params.m_useProxRodRodCollisions = GetBoolOpt( "useProxRodRodCollisions" );
    m_simulation_params.m_useCTRodRodCollisions = GetBoolOpt( "useCTRodRodCollisions" );
    m_simulation_params.m_percentCTRodRodCollisionsAccept = GetScalarOpt( "percentCTRodRodCollisionsAccept" );
    
    m_simulation_params.m_useNonLinearAsFailsafe = GetBoolOpt( "useNonLinearAsFailsafe" );
    m_simulation_params.m_alwaysUseNonLinear = GetBoolOpt( "alwaysUseNonLinear" );
    m_simulation_params.m_maxNewtonIterations = GetIntOpt( "maxNewtonIterations" );

    m_simulation_params.m_simulationManager_limitedMemory = GetBoolOpt( "simulationManager_limitedMemory" );
    m_simulation_params.m_contactProblem_limitedMemory = GetBoolOpt( "contactProblem_limitedMemory" );
    m_simulation_params.m_gaussSeidelTolerance = GetScalarOpt( "gaussSeidelTolerance" );
    m_simulation_params.m_stochasticPruning = GetScalarOpt( "stochasticPruningFraction" );
}

void Scene::setRodCollisionParameters( ElasticStrand& strand )
{
    strand.collisionParameters().m_externalCollisionsRadius = GetScalarOpt( "externalCollisionsRadius" ); //0.05;
    strand.collisionParameters().m_collisionRadius = GetScalarOpt( "collisionRadius" ); //0.0025;
    
    strand.collisionParameters().rootImmunityLength() *= -strand.getTotalRestLength();
    
    strand.collisionParameters().m_frictionCoefficient = GetScalarOpt("strand_mu"); //0.3;
    
    strand.collisionParameters().setAssociatedStrandParameters( strand.getParameters() );
}

void Scene::render( const int& w, const int& h, const int& l, const bool& ct )
{
    clearContacts();

    for( ElasticStrand* sptr = m_strands.begin(); sptr != m_strands.end(); ++sptr )
    {
        m_renderer->render( strand, w, h, l, ct );
    }
    
    for(auto m_itr = m_mesh_renderers.begin(); m_itr != m_mesh_renderers.end(); ++ m_itr)
    {
        (*m_itr)->render();
    }
}

void Scene::getCenter( Vec3d& center )
{
    int render_count = 0;
    for( ElasticStrand* sptr = m_strands.begin(); sptr != m_strands.end(); ++sptr )
    {
        center += m_renderer->calculateObjectCenter( sptr );
        ++render_count;
    }
    
    for(auto m_itr = m_mesh_renderers.begin(); m_itr != m_mesh_renderers.end(); ++ m_itr)
    {
        center += (*m_itr)->calculateObjectCenter();
        ++render_count;
    }
    
    if( render_count > 0 ) center /= render_count;
}

void Scene::getRadius( Scalar& radius, const Vec3d& simCenter )
{
    for( ElasticStrand* sptr = m_strands.begin(); sptr != m_strands.end(); ++sptr )
    {
        const Vec3x center = m_renderer->calculateObjectCenter( sptr );
        const Scalar r = m_renderer->calculateObjectBoundingRadius( sptr );
        radius = std::max( radius, r + ( center - simCenter ).norm() );
    }
    
    for(auto m_itr = m_mesh_renderers.begin(); m_itr != m_mesh_renderers.end(); ++ m_itr)
    {
        const Vec3x center = (*m_itr)->calculateObjectCenter();
        const Scalar r = (*m_itr)->calculateObjectBoundingRadius( center );
        radius = std::max( radius, r + ( center - simCenter ).norm() );
    }
}


// Options :

template <typename T> 
int Scene::AddOption( const std::string& name, const std::string& desc, const T& def );
{
    if( m_options.find(name) != m_options.end() ){
        std::cerr << "Option " << name << " already exists" << std::endl;
        return -1;
    }
    m_options.insert(std::make_pair(name, Option(name, desc, def)));
    return 0;
}

Option* Scene::GetOption(const std::string& name)
{
    if( m_options.find(name) == m_options.end() ){
        std::cerr << "Option " << name << " does not exist\n";
        exit(-1);
    }
    return &(m_options.find(name)->second);
}

bool& Scene::GetBoolOpt(const std::string& name)
{
    return GetOption(name)->b;
}

int& Scene::GetIntOpt(const std::string& name)
{
    return GetOption(name)->i;
}

Scalar& Scene::GetScalarOpt(const std::string& name)
{
    return GetOption(name)->r;
}

Vec3d& Scene::GetVecOpt(const std::string& name)
{
    return GetOption(name)->v;
}

std::string& Scene::GetStringOpt(const std::string& name)
{
    return GetOption(name)->s;
}

int Scene::LoadOptions(const char* filename)
{
    std::ifstream input(filename);
    if (!input.is_open()) {
        std::cerr << "ERROR: File " << filename << " not found\n";
        return -1;
    }
    
    std::string line, option;
    std::istringstream sIn;
    std::string tmp;
    for (getline(input, line); !input.eof(); getline(input, line)) { // FIXME: this won't read the last line of the option file if it is not ended by a newline.
        sIn.clear();
        option.clear();
        sIn.str(line);
        sIn >> option;
        if (option.size() == 0 || option.c_str()[0] == '#') continue;
        std::map<std::string, Option>::iterator itr;
        itr = m_options.find(option);
        if (itr == m_options.end()) {
            std::cerr << "Invalid option: " << option << std::endl;
            continue;
        }
        if (itr->second.type == Option::BOOL) {
            sIn >> tmp;
            if (tmp == "true" || tmp == "1") itr->second.b = true;
            else if (tmp == "false" || tmp == "0") itr->second.b = false;
        } else if (itr->second.type == Option::INT) {
            sIn >> itr->second.i;
        } else if (itr->second.type == Option::SCALAR) {
            sIn >> itr->second.r;
        } else if (itr->second.type == Option::VEC) {
            Vec3d& v = itr->second.v;
            sIn >> v[0];
            sIn >> v[1];
            sIn >> v[2];
        } else if (itr->second.type == Option::STRING) {
            sIn >> itr->second.s;
        } else {
            std::cerr << "Invalid option type" << std::endl;
        }
    }
    input.close();
    
    return 0;
}



void Scene::PrintOptions(std::ostream& os)
{
    std::map<std::string, Option>::const_iterator itr;
    for (itr = m_options.begin(); itr != m_options.end(); ++itr) {
        const Option* opt = &itr->second;
        os << opt->name << " ";
        switch(opt->type) {
            case Option::BOOL:
                if (opt->b) {
                    os << "true";
                } else {
                    os << "false";
                }
                break;
            case Option::INT:
                os << opt->i;
                break;
            case Option::SCALAR:
                os << opt->r;
                break;
            case Option::VEC:
                os << opt->v.format(EIGEN_SPACES_ONLY_IO);
                break;
            case Option::STRING:
                os << opt->s;
                break;
        }
        os << "\t# " << opt->label << std::endl;
    }
}

/// Serialization

void Scene::dumpRods( std::string outputdirectory, int current_frame, int file_width ) const
{
    mkdir(outputdirectory.c_str(), 0755);
    
    std::stringstream name;
    name << std::setfill('0');
    name << outputdirectory << "/rods_" << std::setw(file_width) << current_frame << ".ply";
    std::ofstream os(name.str().c_str());
    
    // header
    os << "ply" << std::endl << "format ascii 1.0" <<std::endl << "comment created by BASim" << std::endl;
    int num_verts = 0;
    for( ElasticStrand* sptr = m_strands.begin(); sptr != m_strands.end(); ++sptr )
    {
        num_verts += sptr->getNumVertices();
    }
    os << "element vertex " << num_verts << std::endl;
    os << "property float x"<< std::endl << "property float y"<< std::endl << "property float z" << std::endl
     << "property int segment" << std::endl << "element face 0"<< std::endl << "property list int int vertex_indices" << std::endl << "end_header " << std::endl;

    // rod vertex positions
    int rod_num = 0;
    for( ElasticStrand* sptr = m_strands.begin(); sptr != m_strands.end(); ++sptr )
    {
        for (int j = 0; j < sptr->getNumVertices(); ++j)
        {
            os << sptr->getVertex(j).x() << " " << sptr->getVertex(j).y() << " " << sptr->getVertex(j).z() << " " << rod_num <<std::endl;
        }
    }
}


template<typename T> void serializeVarHex( const T& var, std::ostream& output_stream )
{
    assert( output_stream.good() );
    T local_var = var;
    output_stream.precision( std::numeric_limits<double>::digits10 + 2);
    output_stream.flags( std::ios_base::fixed | std::ios_base::scientific );
    output_stream.write( reinterpret_cast<char*>( &local_var ), sizeof(T) );
}

template<typename T> void deserializeVarHex( T& var, std::istream& input_stream )
{
    assert( input_stream.good() );
    input_stream.precision( std::numeric_limits<double>::digits10 + 2);
    input_stream.flags( std::ios_base::fixed | std::ios_base::scientific );
    input_stream.read( reinterpret_cast<char*>( &var ), sizeof(T) );
}

void Scene::checkpointSave( std::string outputdirectory ) const
{
    mkdir(outputdirectory.c_str(), 0755);
    
    std::stringstream name;
    name << outputdirectory << "/.chckpnt.bin";
    std::ofstream os( name.str().c_str(), std::ios::out | std::ios::binary | std::ios::trunc );

    // save out sim time
    serializeVarHex( m_t, os );

    os.precision( std::numeric_limits<double>::digits10 + 2);
    os.flags( std::ios_base::fixed | std::ios_base::scientific );

    int numOrRods = 0;

    // rod Xs and dXs
    int rod_num = 0;
    for( ElasticStrand* sptr = m_strands.begin(); sptr != m_strands.end(); ++sptr )
    {

        numOrRods += sptr->getNumVertices() - 1;

        VecXx currDofs = sptr->getCurrentDegreesOfFreedom();
        // if( rod_num == 0 ) std::cout << "currDofs: " << currDofs.transpose() << std::endl;
        VecXx currDisp = sptr->dynamics().getDisplacements();

        for( int v = 0; v < sptr->getNumVertices(); ++v )
        {
            if( v < sptr->getNumVertices() - 1 ){

                serializeVarHex( currDofs[4*v], os );
                serializeVarHex( currDofs[4*v+1], os );
                serializeVarHex( currDofs[4*v+2], os );
                serializeVarHex( currDofs[4*v+3], os );

                serializeVarHex( currDisp[4*v], os );
                serializeVarHex( currDisp[4*v+1], os );
                serializeVarHex( currDisp[4*v+2], os );
                serializeVarHex( currDisp[4*v+3], os );
            }
            else{
                serializeVarHex( currDofs[4*v], os );
                serializeVarHex( currDofs[4*v+1], os );
                serializeVarHex( currDofs[4*v+2], os );

                serializeVarHex( currDisp[4*v], os );
                serializeVarHex( currDisp[4*v+1], os );
                serializeVarHex( currDisp[4*v+2], os );
            }
        }

        Vec3xArray& prevTangCurr =  sptr->getCurrentState().m_referenceFrames1.getPreviousTangents();
        for( int v = 0; v < sptr->getNumVertices() - 1; ++v )
        {
            serializeVarHex( prevTangCurr[v].x(), os );
            serializeVarHex( prevTangCurr[v].y(), os );
            serializeVarHex( prevTangCurr[v].z(), os );
        }

        Vec3xArray& prevTangFut =  sptr->getFutureState().m_referenceFrames1.getPreviousTangents();
        for( int v = 0; v < sptr->getNumVertices() - 1; ++v )
        {
            serializeVarHex( prevTangFut[v].x(), os );
            serializeVarHex( prevTangFut[v].y(), os );
            serializeVarHex( prevTangFut[v].z(), os );
        }

        std::vector<Scalar> currReftwists = sptr->getCurrentReferenceTwistsDirty();
        for( int v = 0; v < sptr->getNumVertices() - 1; ++v )
        {
            serializeVarHex( currReftwists[v], os );
        }

        std::vector<Scalar> FutureReftwists = sptr->getFutureReferenceTwistsDirty();
        for( int v = 0; v < sptr->getNumVertices() - 1; ++v )
        {
            serializeVarHex( FutureReftwists[v], os );
        }

        const Vec3xArray& currRefFrame1 = sptr->getCurrentState().m_referenceFrames1.getDirty();
        for( int v = 0; v < sptr->getNumVertices() - 1; ++v )
        {
            serializeVarHex( currRefFrame1[v].x(), os );
            serializeVarHex( currRefFrame1[v].y(), os );
            serializeVarHex( currRefFrame1[v].z(), os );
        }

        const Vec3xArray& futureRefFrame1 = sptr->getFutureState().m_referenceFrames1.getDirty();
        for( int v = 0; v < sptr->getNumVertices() - 1; ++v )
        {
            serializeVarHex( futureRefFrame1[v].x(), os );
            serializeVarHex( futureRefFrame1[v].y(), os );
            serializeVarHex( futureRefFrame1[v].z(), os );
        }

    }

    serializeVarHex( edgeCounter, os );
    serializeVarHex( numOrRods, os);

    // print out twist bands and all their info
    std::vector< TwistEdge* > tunneledBands = m_strandsManager->m_collisionDetector->m_proxyHistory->tunnelingBands;
    // sort by uniqueID?

    serializeVarHex( tunneledBands.size(), os );

    // std::cout << "edgeCounter: " << edgeCounter << " numOrRods: " << numOrRods << " tunsize: " << tunneledBands.size() << std::endl; 

    for( unsigned t = 0; t < tunneledBands.size(); ++t )
    {
        TwistEdge* edge = tunneledBands[t];

        // os << edge->uniqueID << std::endl;
        serializeVarHex( edge->uniqueID, os );
        // print out parent IDs
        serializeVarHex( edge->parents.first->uniqueID , os );
        serializeVarHex( edge->parents.second->uniqueID, os );

        serializeVarHex( edge->intersectionTwists(), os );
        TwistEdge::printIntersections( edge->intersections, os );
    }

    os.close();
}

static int numOriginalRods;
bool nonOriginal( ElementProxy* ep )
{
    TwistEdge* twist = dynamic_cast< TwistEdge* >( ep );
    //    return twist->uniqueID >= numOriginalRods;
     if( twist ){
         return twist->uniqueID >= numOriginalRods;
     }
     else{
         return false;
     }
    // return (i % 2) != 0;
}

void Scene::checkpointRestore( std::string directory )
{
    GetStringOpt("checkpointDir") = directory;
    std::vector< TwistEdge* >& tunneledBands = m_strandsManager->m_collisionDetector->m_proxyHistory->tunnelingBands;
    tunneledBands.clear();
    // need to clear elementProxies past the original ones....
    checkpointRestore();
}

void Scene::checkpointRestore()
{
    std::string inputdirectory = GetStringOpt("checkpointDir");
    std::stringstream name;
    name << inputdirectory << "/.chckpnt.bin";
    std::ifstream in( name.str().c_str(), std::ios::in | std::ios::binary );
    if( !in.is_open() )
    {
        std::cerr << "ProblemStepper could not find/open checkpoint file: " << name << std::endl;
        exit(1);
    }
    in.precision( std::numeric_limits<double>::digits10 + 2);
    in.flags( std::ios_base::fixed | std::ios_base::scientific );

    // Proper simulation and base rods will already exist

    // need to jump time
    deserializeVarHex( m_t, in );
    m_strandsManager->setTime( m_t );
    std::cout << "m_t: " << m_t << std::endl;
    
    std::cout.precision( std::numeric_limits<double>::digits10 + 2);
    // std::cout.flags( std::ios_base::fixed | std::ios_base::scientific );

    // rod Xs and dXs
    int rod_num = 0;
    for(auto rd_itr = m_rodDatum.begin(); rd_itr != m_rodDatum.end(); ++rd_itr, ++rod_num)
    {
        VecXx DoFs( sptr->getNumVertices() * 4 - 1 );
        VecXx Disps( sptr->getNumVertices() * 4 - 1 );

        // read in each rod completely,
        for( int v = 0; v < sptr->getNumVertices(); ++v )
        {
            if( v < sptr->getNumVertices() - 1 ){
                deserializeVarHex( DoFs[4*v], in );
                deserializeVarHex( DoFs[4*v+1], in );
                deserializeVarHex( DoFs[4*v+2], in );
                deserializeVarHex( DoFs[4*v+3], in );

                deserializeVarHex( Disps[4*v], in );
                deserializeVarHex( Disps[4*v+1], in );
                deserializeVarHex( Disps[4*v+2], in );
                deserializeVarHex( Disps[4*v+3], in );
            }
            else{
                deserializeVarHex( DoFs[4*v], in );
                deserializeVarHex( DoFs[4*v+1], in );
                deserializeVarHex( DoFs[4*v+2], in );

                deserializeVarHex( Disps[4*v], in );
                deserializeVarHex( Disps[4*v+1], in );
                deserializeVarHex( Disps[4*v+2], in );
            }
        }

        Vec3xArray prevTangCurr( sptr->getNumVertices() - 1 );
        for( int v = 0; v < sptr->getNumVertices() - 1; ++v )
        {
            deserializeVarHex( prevTangCurr[v].x(), in );
            deserializeVarHex( prevTangCurr[v].y(), in );
            deserializeVarHex( prevTangCurr[v].z(), in );
        }

        Vec3xArray prevTangFut ( sptr->getNumVertices() - 1 );
        for( int v = 0; v < sptr->getNumVertices() - 1; ++v )
        {
            deserializeVarHex( prevTangFut[v].x(), in );
            deserializeVarHex( prevTangFut[v].y(), in );
            deserializeVarHex( prevTangFut[v].z(), in );
        }

        std::vector<Scalar> currReftwists( sptr->getNumVertices() - 1 ); 
        for( int v = 0; v < sptr->getNumVertices() - 1; ++v )
        {
            deserializeVarHex( currReftwists[v], in );
        }

        std::vector<Scalar> FutureReftwists( sptr->getNumVertices() - 1 ); 
        for( int v = 0; v < sptr->getNumVertices() - 1; ++v )
        {
            deserializeVarHex( FutureReftwists[v], in );
        }

        Vec3xArray currRefFrame1( sptr->getNumVertices() - 1 );
        for( int v = 0; v < sptr->getNumVertices() - 1; ++v )
        {
            deserializeVarHex( currRefFrame1[v].x(), in );
            deserializeVarHex( currRefFrame1[v].y(), in );
            deserializeVarHex( currRefFrame1[v].z(), in );
        }

        Vec3xArray futureRefFrame1( sptr->getNumVertices() - 1 );
        for( int v = 0; v < sptr->getNumVertices() - 1; ++v )
        {
            deserializeVarHex( futureRefFrame1[v].x(), in );
            deserializeVarHex( futureRefFrame1[v].y(), in );
            deserializeVarHex( futureRefFrame1[v].z(), in );
        }

        sptr->dynamics().setDisplacements( Disps );
        sptr->setCurrentDegreesOfFreedom( DoFs );
        sptr->setFutureDegreesOfFreedom( DoFs );
        sptr->getCurrentState().m_referenceFrames1.setPreviousTangents( prevTangCurr );
        sptr->getFutureState().m_referenceFrames1.setPreviousTangents( prevTangFut );
        sptr->getCurrentState().m_referenceFrames1.cleanSet( currRefFrame1 );
        sptr->getFutureState().m_referenceFrames1.cleanSet( futureRefFrame1 );
        sptr->setCurrentReferenceTwistsClean( currReftwists );
        sptr->setFutureReferenceTwistsClean( FutureReftwists );
    }

    deserializeVarHex( edgeCounter, in );
    int numOrRods;
    deserializeVarHex( numOrRods, in);

    unsigned numTunneledBands;
    deserializeVarHex( numTunneledBands, in );

    numOriginalRods = numOrRods;

    std::vector< TwistEdge* >& tunneledBands = m_strandsManager->m_collisionDetector->m_proxyHistory->tunnelingBands;
    std::vector<ElementProxy*>& originalTE = m_strandsManager->m_collisionDetector->m_elementProxies;

    // clear existing tunnelingBands, and clear all elementproxies past original by index (so when we save->restore we dont create duplicates)
    tunneledBands.clear();
    for( unsigned w = 0; w < originalTE.size(); ++w ){
        TwistEdge* twist = dynamic_cast< TwistEdge* >( originalTE[w] );
        if( twist && !twist->isTwistedBand ){
            twist->children.clear();
        }
    }   

    // need to clear elementproxies past the original ones...
    originalTE.erase( std::remove_if(originalTE.begin(), originalTE.end(), nonOriginal), originalTE.end() );

    deserializeVarHex( numOrRods, in );

    for( unsigned t = 0; t < numTunneledBands; ++t )
    {
        int uID;    
        // read in parent uniqueIDs and then look for them in the regular
        int alphaID, betaID;

        deserializeVarHex( uID, in );
        deserializeVarHex( alphaID, in );
        deserializeVarHex( betaID, in );

        TwistEdge* alpha = NULL;
        TwistEdge* beta = NULL;

        for( unsigned e = 0; e < originalTE.size(); ++e )
        {
            TwistEdge* twist = dynamic_cast< TwistEdge* >( originalTE[e] );
            if( !twist ) continue;
            if( twist->uniqueID == alphaID ) alpha = twist;
            if( twist->uniqueID == betaID ) beta = twist;
            if( alpha != NULL && beta != NULL ) break;
        }

        if( alpha == NULL || beta == NULL ){
            std::cerr << "could not find parent edges, exiting. " << std::endl; std::exit(EXIT_FAILURE);
        }

        TwistEdge* te = new TwistEdge( alpha, beta, uID );

        // after reading in parents, use them to update children with new TwistEdge
        alpha->children.push_back( std::pair< int, TwistEdge* >( beta->uniqueID, te ) );
        beta->children.push_back( std::pair< int, TwistEdge* >( alpha->uniqueID, te ) );

        int intersectionTwists;
        deserializeVarHex( intersectionTwists, in );

        for( int it = 0; it < intersectionTwists; ++it ){
            Scalar oAngle;
            Scalar cAngle;
            int cTwists;

            deserializeVarHex( oAngle, in );
            deserializeVarHex( cAngle, in );
            deserializeVarHex( cTwists, in );

            TwistIntersection* newIntr = new TwistIntersection( oAngle );
            newIntr->currAngle = cAngle;
            newIntr->coplanarTwists = cTwists;
            te->intersections.push( newIntr );
        }

        originalTE.push_back( te );
        tunneledBands.push_back( te );
    } 

    in.close();
}


