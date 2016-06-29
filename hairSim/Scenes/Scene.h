#ifndef SCENE_H
#define SCENE_H

#include "../../StrandSim/Core/Definitions.hh"
#include "../../StrandSim/Collision/CollisionUtils.hh"
#include "../../StrandSim/Dynamic/StrandDynamicTraits.hh"
#include "../../StrandSim/Dynamic/DynamicStepper.hh"
#include "../../StrandSim/Dynamic/DOFScriptingController.hh"
#include "../../StrandSim/Utils/Distances.hh"
#include "../../StrandSim/Utils/Timer.hh"

#include "Control/SimpleMeshController.hh"

#include "Core/ElasticStrand.h"
#include "Dynamic/StrandManager.h"
#include "Render/MeshRenderer.h"
#include "Scenes/Option.h"

class Scene
{
public:
    
    explicit Scene( const std::string& name = "", const std::string& desc = "" );
    virtual ~Scene();
    
    void setup();
    bool step();
    void render( const int& w, const int& h, const int& l, const bool& ct );

    // Scene Options:
        
    void PrintOptions(std::ostream& os);

    template <typename T>
    int AddOption( const std::string& name, const std::string& desc, const T& def );
   
    int LoadOptions( const char* filename );
    int LoadOptions( const std::string& filename ){ return LoadOptions( filename.c_str() ); }

    Option* GetOption( const std::string& name );
    bool& GetBoolOpt( const std::string& name );
    int& GetIntOpt( const std::string& name );
    Scalar& GetScalarOpt( const std::string& name );
    Vec3d& GetVecOpt( const std::string& name );
    std::string& GetStringOpt( const std::string& name );

    // Util functions:

    void isSimulated( bool isSimulated ){ m_isSimulated = isSimulated; }
    bool isSimulated() const { return m_isSimulated; }

    void getCenter( Vec3d& center );
    void getRadius( Scalar& radius, const Vec3d& center );
    
    void setSimulationParameters();
        
    Scalar getTime(){ return m_t; }
    void setTime( Scalar t){ m_t = t; }
    
    Scalar getDt(){ return m_dt; }
    void setDt( Scalar dt) { m_dt = dt; }
        
    // output
    void dumpRods( std::string outputdirectory, int current_frame, int file_width ) const;
    void checkpointSave( std::string outputdirectory ) const;
    void checkpointRestore( std::string directory );
    void checkpointRestore();

    std::string m_problemName;
    std::string m_problemDesc;    
    
protected:

    friend class SceneUtils;
    
    virtual void setupStrands() = 0;
    virtual void setupMeshes() = 0;    
    virtual bool executeScript() = 0; // execute scripting and any per-step updates

    void setRodCollisionParameters( ElasticStrand& strand );
    void addOptions();
    void clearContacts();

    ///////////////////////////

    Scalar m_t;
    Scalar m_dt;
    bool m_isSimulated;

    SimulationParameters m_simulation_params;
    Simulation* m_strandsManager;
    std::map< std::string, Option > m_options;
    std::vector< ElasticStrand* > m_strands;
    std::vector< TriMesh* > m_meshes;
    StrandRenderer* m_renderer;

};
#endif
