#ifndef SCENE_UTILS_H
#define SCENE_UTILS_H

class SceneUtils
{
public:
      
    // hair generation
    static void findOrthogonal( Vec3x& v, const Vec3x& u );
    static void genCurlyHair( const Vec3x& initnorm, const Vec3x& startpoint, const double& dL, 
        const int& nv, std::vector<Vec3x>& vertices, double curl_radius, double curl_density, double root_length );
    static void genStraightHair( const Vec3x& initnorm, const Vec3x& startpoint, const double& dL,
        const int& nv, std::vector<Vec3d>& vertices, double root_length );
    
    // transform utilities
    static void transformTriangleObject( TriangularMesh& triMesh, Mat3x& transformation, Vec3x& center, Vec3x& translate );
    static void freezeTriangleObject( TriangularMesh& triMesh );
    static void transformRodRoot( ElasticStrand* strand, Mat3x& transformation, Vec3x& center, Vec3x& translate, int substeps );
    static void transformRodRootVtx( ElasticStrand* strand, Mat3x& transformation, Vec3x& center, Vec3x& translate, int vtx, int substeps );
    static void translateRodVertex( ElasticStrand* strand, int vtx_id, Vec3x& translate, int substeps );
    static void freezeRodRoot( ElasticStrand* strand );
    static void freezeVertex( ElasticStrand* strand, int vtx );

    // output
    static void dumpMesh( std::string outputdirectory, int current_frame, int file_width ) const;

protected:

private:

};

#endif
