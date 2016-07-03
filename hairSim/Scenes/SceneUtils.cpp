#include "SceneUtils.h"
#include "Scene.h"
#include "../Utils/Option.h"
#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"
#include <sys/stat.h>
#include <iomanip>
#include <fstream>
#define EPSILON 1.0e-12

void SceneUtils::findOrthogonal( Vec3& v, const Vec3& u )
{
    assert(u.norm() != 0);
    
    v.setZero();
    int max = 0;
    for (int i = 0; i < u.size(); ++i) {
        if (u[i] == 0) {
            v[i] = 1;
            return;
        }
        if (fabs(u[i]) > fabs(u[max])) max = i;
    }
    
    int idx = (max + 1) % u.size();
    v[idx] = u[max];
    v[max] = -u[idx];
    v.normalize();
    
    assert(std::abs(u.dot(v)) < EPSILON );
}

void SceneUtils::genCurlyHair( const Vec3& initnorm, const Vec3& startpoint, const double& dL, const int& nv, std::vector<Vec3>& vertices, double curl_radius, double curl_density, double root_length, double curl_density_perturbation, double curl_radius_perturbation )
{
    // generate an orthonormal frame
    Vec3 p1;
    findOrthogonal( p1, initnorm );
    Vec3 p2 = p1.cross( initnorm );
    
    Scalar curl_density_eps = curl_density_perturbation;
    Scalar curl_radius_eps = curl_radius_perturbation;
  
    if ( curl_density_eps > 0 && curl_radius_eps > 0)
    {
        typedef boost::minstd_rand rng_type;
        typedef boost::uniform_real<double> dist_type;
        typedef boost::variate_generator<rng_type&, dist_type> generator_type;
        dist_type dist_density( 1. , 1. + curl_density_eps );
        dist_type dist_radius( 1. , 1. + curl_radius_eps );
        
        struct timeval cur ;
        gettimeofday(&cur, 0) ;
        rng_type generator(cur.tv_usec);
        
        generator_type gen_density(generator, dist_density);
        generator_type gen_radius(generator, dist_radius);
        
        curl_density *= gen_density();
        curl_radius *= gen_radius();
    }
  
    vertices.push_back( startpoint );
    Vec3 freepoint( startpoint + root_length * initnorm );
    vertices.push_back( freepoint );
    // for each curve sample vertex
    Scalar xa = M_PI/( curl_density * 4 ); // 0 // start curve parameter
    Scalar xb = 0; // end curve parameter
    for( int j = 1; j < nv-1; ++j )
    {
        xb = ( dL + xa + curl_radius * curl_radius * curl_density * curl_density * xa)/(1 + curl_radius * curl_radius * curl_density * curl_density ); //upate to get length dL along curve
        vertices.push_back(
                           freepoint
                           + xb * initnorm
                           +  curl_radius * cos( xb * curl_density  ) * p2
                           +  curl_radius * sin( xb * curl_density  ) * p1
                           );
        xa = xb; // next...
    }
}

void SceneUtils::genStraightHair( const Vec3& initnorm, const Vec3& startpoint, const double& dL, const int& nv, std::vector<Vec3>& vertices, double root_length  )
{
    vertices.push_back( startpoint );
    vertices.push_back( startpoint + root_length * initnorm );
    for( int j = 0; j < nv-2; ++j ){
        vertices.push_back( vertices.back() + dL * initnorm );
    }
}

void SceneUtils::transformTriangleObject( TriMesh& triMesh, Mat3x& transformation, Vec3& center, Vec3& translate )
{
    for( unsigned i = 0; i < triMesh.nv(); ++i )
    {
        Vec3 vert = triMesh.getVertex(i);
        Vec3 vertNext = transformation * (vert - center) + center + translate;
        triMesh.setVertex(i, vertNext);
        triMesh.setDisplacement( i, (vertNext-vert) );
    }
}

void SceneUtils::freezeTriangleObject( TriMesh& triMesh )
{
    Vec3 disp( 0. , 0. , 0. ) ;
    for( unsigned i = 0 ; i < triMesh.nv() ; ++i )
    {
        triMesh.setDisplacement( i, disp );
    }
}

void SceneUtils::transformRodRoot( ElasticStrand* strand, Mat3x& transformation, Vec3& center, Vec3& translate )
{
    for(size_t vtx = 0; vtx < 2; ++vtx )
    {
        Vec3 vert = strand->getVertex( vtx );
        Vec3 vertNext = transformation * (vert - center) + center + translate;
        strand->dynamics().getScriptingController()->setVertexDisplacement( vtx, (vertNext - vert) );
    }
}

void SceneUtils::transformRodRootVtx( ElasticStrand* strand, Mat3x& transformation, Vec3& center, Vec3& translate, int vtx )
{
    Vec3 vert = strand->getVertex( vtx );
    Vec3 vertNext = transformation * (vert - center) + center + translate;
    strand->dynamics().getScriptingController()->setVertexDisplacement( vtx, (vertNext - vert) );
    strand->dynamics().getScriptingController()->setThetaDisplacement( 0, 0. );
}

void SceneUtils::translateRodVertex( ElasticStrand* strand, int vtx_id, Vec3& translate )
{
    strand->dynamics().getScriptingController()->setVertexDisplacement( vtx_id, translate );
}

void SceneUtils::freezeRodRoot( ElasticStrand* strand )
{
    strand->dynamics().getScriptingController()->freezeVertices( 0, true );
    strand->dynamics().getScriptingController()->freezeVertices( 1 );
}

void SceneUtils::freezeVertex( ElasticStrand* strand, int vtx )
{
    strand->dynamics().getScriptingController()->freezeVertices( vtx );
}

void SceneUtils::dumpMesh( std::string outputdirectory, int current_frame, int file_width, const std::vector< TriMesh* >& meshes )
{
    mkdir(outputdirectory.c_str(), 0755);
    int mesh_num = 0;
    for( auto m_itr = meshes.begin(); m_itr != meshes.end(); ++ m_itr, ++ mesh_num )
    {
        // new obj per mesh
        std::stringstream name;
        name << std::setfill('0');
        name << outputdirectory << "/mesh" <<  mesh_num << "_" << std::setw(file_width) << current_frame << ".obj";
        std::ofstream os(name.str().c_str());
        // header
        os << "# obj created by HairSim" << std::endl;
        // vertices
        TriMesh* mesh = *m_itr;
        for ( size_t v = 0; v < mesh->nv(); ++v )
        {
            os << "v " << mesh->getVertex(v).x() << " " << mesh->getVertex(v).y() << " " << mesh->getVertex(v).z() << std::endl;
        }
        // triangles
        for(auto f_itr = mesh->getFaces().begin(); f_itr != mesh->getFaces().end(); ++ f_itr)
        {
            os << "f " << ( *f_itr ).idx[0] + 1 << " " << ( *f_itr ).idx[1] + 1 << " " << ( *f_itr ).idx[2] + 1 << std::endl;
        }
        os.close();
    }
}

