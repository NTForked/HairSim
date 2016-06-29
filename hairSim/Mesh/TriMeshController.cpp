#include "TriangularMesh.h"

// Duration of guessed normal sign validity, in number of frames without query
#define NORMAL_SIGN_VALIDITY 5

using namespace std;

SimpleMeshController::SimpleMeshController( double i_time, double i_dt ):
    MeshScriptingController( i_time, i_dt ), 
    m_isStaticMesh( true ),
    m_startMeshTime( 0. ), 
    m_endMeshTime( 0. ), 
    m_lastExecutionTime( 0 ), 
    m_defaultFrictionCoefficient( 0. ), 
{
    m_startTime = i_time;
    m_mesh = new TriMesh();
    m_mesh->setAssociatedController( this );
}

SimpleMeshController::~SimpleMeshController()
{
    delete m_mesh;
    m_mesh = NULL;
}

void SimpleMeshController::isStaticMesh( const bool i_isStaticMesh )
{
    m_isStaticMesh = i_isStaticMesh;
}

bool SimpleMeshController::execute( bool updateLevelSet )
{
    return true;
}

short SimpleMeshController::knowsNormalSign( bool atPreviousStep, unsigned faceIndex,
                                            unsigned rodIndex, unsigned vertex )
{
    return 1;
}

bool SimpleMeshController::loadMesh( std::string& obj_file_name )
{
    // intialize all triangular meshes
    ObjParser objparser;
    objparser.loadTriangularMesh( obj_file_name, *m_mesh);
    
    std::cout<< "# Loaded mesh: " << obj_file_name << std::endl;
    return true; 
}
