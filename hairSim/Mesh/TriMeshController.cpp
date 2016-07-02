#include "TriMesh.h"
#include "TriMeshController.h"

// Duration of guessed normal sign validity, in number of frames without query
#define NORMAL_SIGN_VALIDITY 5

using namespace std;

TriMeshController::TriMeshController( double i_time, double i_dt ):
    m_time( i_time),
    m_dt( i_dt ), 
    m_isStaticMesh( true ),
    m_startMeshTime( 0. ), 
    m_endMeshTime( 0. ), 
    m_lastExecutionTime( 0 ), 
    m_defaultFrictionCoefficient( 0. )
{
    m_startTime = i_time;
    m_mesh = new TriMesh();
    m_mesh->setController( this );
}

TriMeshController::~TriMeshController()
{
    delete m_mesh;
    m_mesh = NULL;
}

void TriMeshController::isStaticMesh( const bool i_isStaticMesh )
{
    m_isStaticMesh = i_isStaticMesh;
}

bool TriMeshController::execute( bool updateLevelSet )
{
    return true;
}

short TriMeshController::knowsNormalSign( bool atPreviousStep, unsigned faceIndex,
                                            unsigned rodIndex, unsigned vertex )
{
    return 1;
}

bool TriMeshController::loadMesh( std::string& obj_file_name )
{
    // intialize all triangular meshes
    ObjParser objparser;
    objparser.loadTriMesh( obj_file_name, *m_mesh);
    
    std::cout<< "# Loaded mesh: " << obj_file_name << std::endl;
    return true; 
}
