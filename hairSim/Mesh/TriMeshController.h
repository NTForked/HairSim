#ifndef TriMeshController_HH_
#define SIMPLESHCONTROLLER_HH_

#include "../Utils/Definitions.h"
#include "ObjParser.h"

#include <map>


/* H
    should be changed to a struct/ simple class
    which holds some relevant info/params/utils for 
    triMesh class

    maybe put the Scene tri-mesh utils here?
    (scripting positions etc)
*/

class TriMesh;

class TriMeshController
{
public:
    TriMeshController( double i_time, double i_dt );

    virtual ~TriMeshController();
    
    virtual void isStaticMesh( const bool i_isStaticMesh );

    virtual bool execute( bool updateLevelSet );

    const TriMesh* getMesh() const
    {
        return m_mesh;
    }

    double getDefaultFrictionCoefficient() const
    {
        return m_defaultFrictionCoefficient;
    }

    void setDefaultFrictionCoefficient( const double frictionCoeff )
    {
        m_defaultFrictionCoefficient = frictionCoeff;
    }

    virtual const std::vector<double>& getFrictionCoefficients() const
    {
        return m_frictionCoefficients;
    }

    virtual std::vector<double>& getFrictionCoefficients()
    {
        return m_frictionCoefficients;
    }
    
    bool loadMesh( std::string& obj_file_name );

    bool hasLevelSet() const
    {
        return false;
    }
    
    const std::vector<bool>& getEnabledVertices() const
    {
        return m_enabledVertices;
    }
    
    virtual std::vector<bool>& getEnabledVertices()
    {
        return m_enabledVertices;
    }
    
    double getLevelSetForceThickness() const
    {
        std::cout << "Level sets not supported by TriMeshController\n";
        exit(1);
    }
    
    double getLevelSetForceStrength() const
    {
        std::cout << "Level sets not supported by TriMeshController\n";
        exit(1);
    }
    
    TriMesh* getMesh()
    {
        return m_mesh;
    }

   // Virtual functions that should be inherited for meshes that are collidable on both sides
    virtual bool collideOnBothSides() const
    {
        return false;
    }

    // Whether the controller is sure of the side of which a vertex should stay. Returns:
    // param atPreviousStep ; if true, do not tka einto account guesses from the current step
    //  1 : normal is known and correspond to counter-clockwise vector product of vertices
    //  0 : normal is unknown
    // -1 : normal is known and correspond to clockwise vector product of vertices
    virtual short knowsNormalSign( bool atPreviousStep, unsigned faceIndex, unsigned rodIndex,
            unsigned vertex );

    virtual void setNormalSign( short , float , unsigned , unsigned ,
            unsigned  )
    {}


private:
    
    Scalar m_time;
    Scalar m_dt;

    bool m_isStaticMesh;
    TriMesh* m_mesh;

    double m_startMeshTime;
    double m_endMeshTime;

    double m_startTime;
    double m_lastExecutionTime;

    std::vector<bool> m_enabledVertices;
    std::vector<double> m_frictionCoefficients;
    double m_defaultFrictionCoefficient;
};

#endif
