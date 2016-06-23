#include "SceneUtils.h"

void SceneUtils::findOrthogonal( Vec3x& v, const Vec3x& u )
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

void SceneUtils::generateCurlyHair( const Vec3x& initnorm, const Vec3x& startpoint, const double& dL, const int& nv, std::vector<Vec3x>& vertices, double curl_radius, double curl_density, double root_length )
{
    // generate an orthonormal frame
    Vec3x p1;
    findOrthogonal( p1, initnorm );
    Vec3x p2 = p1.cross( initnorm );
    
    Scalar curl_density_eps = GetScalarOpt( "curl_density_perturbation" );
    Scalar curl_radius_eps = GetScalarOpt( "curl_radius_perturbation" );
  
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
    Vec3x freepoint( startpoint + root_length * initnorm );
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

void SceneUtils::generateStraightHair( const Vec3x& initnorm, const Vec3x& startpoint, const double& dL, const int& nv, std::vector<Vec3x>& vertices, double root_length  )
{
    vertices.push_back( startpoint );
    vertices.push_back( startpoint + root_length * initnorm );
    for( int j = 0; j < nv-2; ++j ){
        vertices.push_back( vertices.back() + dL * initnorm );
    }
}

void SceneUtils::transformTriangleObject( TriangularMesh& triMesh, Mat3x& transformation, Vec3x& center, Vec3x& translate )
{
    for( unsigned i = 0; i < triMesh.nv(); ++i )
    {
        Vec3x vert = triMesh.getVertex(i);
        Vec3x vertNext = transformation * (vert - center) + center + translate;
        triMesh.setVertex(i, vertNext);
        triMesh.setDisplacement( i, (vertNext-vert) );
    }
}

void SceneUtils::freezeTriangleObject( TriangularMesh& triMesh )
{
    Vec3x disp( 0. , 0. , 0. ) ;
    for( unsigned i = 0 ; i < triMesh.nv() ; ++i )
    {
        triMesh.setDisplacement( i, disp );
    }
}

void SceneUtils::transformRodRoot( ElasticStrand* strand, Mat3x& transformation, Vec3x& center, Vec3x& translate )
{
    for(size_t vtx = 0; vtx < 2; ++vtx )
    {
        Vec3x vert = strand->getVertex( vtx );
        Vec3x vertNext = transformation * (vert - center) + center + translate;
        strand->dynamics()->getDofController().setVertexDisplacement( vtx, (vertNext - vert) / m_subSteps );
    }
}

void SceneUtils::transformRodRootVtx( ElasticStrand* strand, Mat3x& transformation, Vec3x& center, Vec3x& translate, int vtx )
{
    Vec3x vert = strand->getVertex( vtx );
    Vec3x vertNext = transformation * (vert - center) + center + translate;
    strand->dynamics()->getDofController().setVertexDisplacement( vtx, (vertNext - vert) / m_subSteps );
    strand->dynamics()->getDofController().setThetaDisplacement( 0, 0. );
}

void SceneUtils::translateRodVertex( ElasticStrand* strand, int vtx_id, Vec3x& translate )
{
    strand->dynamics()->getDofController().setVertexDisplacement( vtx_id, translate / m_subSteps );
}

void SceneUtils::freezeRodRoot( ElasticStrand* strand )
{
    strand->dynamics()->getDofController().freezeRootVertices<2>();
}

void SceneUtils::freezeVertex( ElasticStrand* strand, int vtx )
{
    strand->dynamics()->getDofController().freezeVertices( vtx );
}

void SceneUtils::dumpMesh( std::string outputdirectory, int current_frame, int file_width ) const
{
    mkdir(outputdirectory.c_str(), 0755);
    int mesh_num = 0;
    for(auto m_itr = m_meshScripting_controllers.begin(); m_itr != m_meshScripting_controllers.end(); ++ m_itr, ++ mesh_num)
    {
        // new obj per mesh
        std::stringstream name;
        name << std::setfill('0');
        name << outputdirectory << "/mesh" <<  mesh_num << "_" << std::setw(file_width) << current_frame << ".obj";
        std::ofstream os(name.str().c_str());
        // header
        os << "# obj created by StrandSim" << std::endl;
        // vertices
        strandsim::TriangularMesh* mesh = (*m_itr)->getCurrentMesh();
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

