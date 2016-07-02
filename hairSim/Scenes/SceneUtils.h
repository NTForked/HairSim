#ifndef SCENE_UTILS_H
#define SCENE_UTILS_H

#include "../Utils/Definitions.h"
#include "../Mesh/TriMesh.h"
#include "../Strand/ElasticStrand.h"

class SceneUtils
{
public:
      
    // hair generation
    static void findOrthogonal( Vec3& v, const Vec3& u );
    static void genCurlyHair( const Vec3& initnorm, const Vec3& startpoint, const double& dL, 
        const int& nv, std::vector<Vec3>& vertices, double curl_radius, double curl_density, double root_length, double curl_density_perturbation, double curl_radius_perturbation  );
    static void genStraightHair( const Vec3& initnorm, const Vec3& startpoint, const double& dL,
        const int& nv, std::vector<Vec3d>& vertices, double root_length );
    
    // transform utilities
    static void transformTriangleObject( TriMesh& triMesh, Mat3x& transformation, Vec3& center, Vec3& translate );
    static void freezeTriangleObject( TriMesh& triMesh );
    static void transformRodRoot( ElasticStrand* strand, Mat3x& transformation, Vec3& center, Vec3& translate );
    static void transformRodRootVtx( ElasticStrand* strand, Mat3x& transformation, Vec3& center, Vec3& translate, int vtx );
    static void translateRodVertex( ElasticStrand* strand, int vtx_id, Vec3& translate );
    static void freezeRodRoot( ElasticStrand* strand );
    static void freezeVertex( ElasticStrand* strand, int vtx );

    // output
    static void dumpMesh( std::string outputdirectory, int current_frame, int file_width, const std::vector< TriMesh* >& meshes );

protected:

private:

};

#endif
