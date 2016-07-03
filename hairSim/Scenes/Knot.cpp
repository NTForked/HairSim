#include "Knot.h"

#define PI 3.14159265358979323846

Knot::Knot() :
Scene("Knot", "TrefoilKnot, overhand knot pulled tight"),
m_radius(3.)
{
    AddOption( m_problemName, m_problemDesc, "" );

    // Global opts
    GetScalarOpt("dt") = 1e-4;
    GetVecOpt("gravity") = Vec3( 0.0, 0.0, 0.0 );
    AddOption("translation","how much to move the nonfixed end per timestep", Vec3( 0.0, -0.0001, 0.0 ) ); //   0.25, 0, 0.5
    
    // Rod opts
    GetIntOpt("nv") = 84;//35 + 20; //try 11, 981.0, uneven vertices (shows it better)
    
    // Additional parameters:
    AddOption("totalLength","enforced strand length", 35.0 );
    //  AddOption("thickness", "thickness", 0.5);
    //  AddOption("stiffness", "stiffness", 1000.0);
    //  AddOption("mass-damping", "mass damping for the rod", 0.0);

    // Pre-setup to default values:
    GetScalarOpt( "stochasticPruningFraction" ) = 0.5;
    GetBoolOpt("useProxRodRodCollisions") = true; // skip Proximity == true.
    GetScalarOpt("collisionRadius") = 0.005;
    GetBoolOpt("useCTRodRodCollisions") = true;
    
    GetScalarOpt("strand_mu") = 0.0;

}

Knot::~Knot()
{}

void Knot::setupStrands()
{
    std::cout << "Scene:: SingleContact" << std::endl;

    // discrete rod params
    const int nVertices = GetIntOpt("nv");
    const int subdivision = 5;

    // Prepare initial rod/strand position
    Vec3Array i_vertices;


    const int nDOFs = 4 * nVertices - 1;
    
    // rod params
    const Scalar totalLength = GetScalarOpt("totalLength");
    
    
//    const Scalar radiusA = 0.000254;
//    const Scalar youngsModulus = 67.5e9;
//    const Scalar shearModulus = 22.5e9;
//    const Scalar density = 6450;
    
    const Scalar radiusA = GetScalarOpt("radius");
    const Scalar youngsModulus = GetScalarOpt("youngs-modulus");
    const Scalar shearModulus = GetScalarOpt("shear-modulus");
    const Scalar density = GetScalarOpt("density");
    
    const Scalar airDrag = GetScalarOpt("air-drag");
    const Scalar viscosity = GetScalarOpt("viscosity");
    const Scalar baseRotation = GetScalarOpt("base-rotation");

//    i_vertices.push_back( Vec3d( -0.560566, 1, 7.1418));
//    i_vertices.push_back( Vec3d( -0.540374, 1, 6.9145));
//    i_vertices.push_back( Vec3d( -0.51952, 1, 6.68592));
//    i_vertices.push_back( Vec3d( -0.4976, 1, 6.45506));
//    i_vertices.push_back( Vec3d( -0.474659, 1, 6.22119));
//    i_vertices.push_back( Vec3d( -0.451276, 1, 5.9));
    i_vertices.push_back( Vec3d( -0.40653, 0.1, 5.49881));
    i_vertices.push_back( Vec3d( -0.37669, 0.1, 4.9966));
    i_vertices.push_back( Vec3d( -0.368871, 0.1, 4.47473));
    i_vertices.push_back( Vec3d( -0.379639, 0.1, 3.93111));
    i_vertices.push_back( Vec3d( -0.400419, 0.0714052, 3.3684));
    i_vertices.push_back( Vec3d( -0.415871, 0.0483848, 2.91412));
    i_vertices.push_back( Vec3d( -0.420971, 0.0812021, 2.46022));
    i_vertices.push_back( Vec3d( -0.412144, 0.163211, 2.01264));
    i_vertices.push_back( Vec3d( -0.391751, 0.278497, 1.57283));
    i_vertices.push_back( Vec3d( -0.368333, 0.404173, 1.13604));
    i_vertices.push_back( Vec3d( -0.353446, 0.531234, 0.622013));
    i_vertices.push_back( Vec3d( -0.347984, 0.602504, 0.0971484));
    i_vertices.push_back( Vec3d( -0.32824, 0.606284, -0.432171));
    i_vertices.push_back( Vec3d( -0.283765, 0.538083, -0.955573));
    i_vertices.push_back( Vec3d( -0.215297, 0.392705, -1.4603));
    i_vertices.push_back( Vec3d( -0.177059, 0.29022, -1.69638));
    i_vertices.push_back( Vec3d( -0.137769, 0.16692, -1.9221));
    i_vertices.push_back( Vec3d( -0.0995085, 0.0227982, -2.13531));
    i_vertices.push_back( Vec3d( -0.0643051, -0.141623, -2.33386));
    i_vertices.push_back( Vec3d( -0.033965, -0.325247, -2.51568));
    i_vertices.push_back( Vec3d( -0.00842831, -0.541767, -2.68999));
    i_vertices.push_back( Vec3d( 0.0086625, -0.776033, -2.84079));
    i_vertices.push_back( Vec3d( 0.0167866, -1.02519, -2.96636));
    i_vertices.push_back( Vec3d( 0.0161079, -1.28619, -3.06532));
    i_vertices.push_back( Vec3d( 0.00743748, -1.55592, -3.13667));
    i_vertices.push_back( Vec3d( -0.0139794, -1.91836, -3.18746));
    i_vertices.push_back( Vec3d( -0.0425906, -2.28386, -3.18942));
    i_vertices.push_back( Vec3d( -0.0736473, -2.64627, -3.14362));
    i_vertices.push_back( Vec3d( -0.10262, -3.00019, -3.05243));
    i_vertices.push_back( Vec3d( -0.125862, -3.34101, -2.91935));
    i_vertices.push_back( Vec3d( -0.140973, -3.66402, -2.74926));
    i_vertices.push_back( Vec3d( -0.147017, -3.96771, -2.54618));
    i_vertices.push_back( Vec3d( -0.144036, -4.25013, -2.31438));
    i_vertices.push_back( Vec3d( -0.132742, -4.50981, -2.05761));
    i_vertices.push_back( Vec3d( -0.114149, -4.74555, -1.77907));
    i_vertices.push_back( Vec3d( -0.0891148, -4.95784, -1.47879));
    i_vertices.push_back( Vec3d( -0.05898, -5.14322, -1.16165));
    i_vertices.push_back( Vec3d( -0.0251582, -5.30056, -0.830041));
    i_vertices.push_back( Vec3d( 0.010555, -5.42913, -0.486454));
    i_vertices.push_back( Vec3d( 0.045903, -5.52887, -0.13338));
    i_vertices.push_back( Vec3d( 0.0858806, -5.61407, 0.318084));
    i_vertices.push_back( Vec3d( 0.115689, -5.65832, 0.776158));
    i_vertices.push_back( Vec3d( 0.130615, -5.66525, 1.23703));
    i_vertices.push_back( Vec3d( 0.127749, -5.63764, 1.69738));
    i_vertices.push_back( Vec3d( 0.10661, -5.57601, 2.15391));
    i_vertices.push_back( Vec3d( 0.0710824, -5.48271, 2.58415));
    i_vertices.push_back( Vec3d( 0.0240035, -5.35127, 3.00317));
    i_vertices.push_back( Vec3d( -0.0302718, -5.1757, 3.4048));
    i_vertices.push_back( Vec3d( -0.0872761, -4.95037, 3.78035));
    i_vertices.push_back( Vec3d( -0.143345, -4.67204, 4.11868));
    i_vertices.push_back( Vec3d( -0.195518, -4.34851, 4.40205));
    i_vertices.push_back( Vec3d( -0.245034, -3.98089, 4.62586));
    i_vertices.push_back( Vec3d( -0.29355, -3.57944, 4.78132));
    i_vertices.push_back( Vec3d( -0.342537, -3.15702, 4.86413));
    i_vertices.push_back( Vec3d( -0.390972, -2.72664, 4.87496));
    i_vertices.push_back( Vec3d( -0.422388, -2.41994, 4.84121));
    i_vertices.push_back( Vec3d( -0.448459, -2.11804, 4.77509));
    i_vertices.push_back( Vec3d( -0.466561, -1.82374, 4.6789));
    i_vertices.push_back( Vec3d( -0.474799, -1.53965, 4.55474));
    i_vertices.push_back( Vec3d( -0.472653, -1.26844, 4.40429));
    i_vertices.push_back( Vec3d( -0.462215, -1.02785, 4.24021));
    i_vertices.push_back( Vec3d( -0.445772, -0.803347, 4.05515));
    i_vertices.push_back( Vec3d( -0.426239, -0.597026, 3.8503));
    i_vertices.push_back( Vec3d( -0.406725, -0.410587, 3.62719));
    i_vertices.push_back( Vec3d( -0.389998, -0.245201, 3.38785));
    i_vertices.push_back( Vec3d( -0.376269, -0.0730787, 3.07809));
    i_vertices.push_back( Vec3d( -0.372348, 0.0674958, 2.75252));
    i_vertices.push_back( Vec3d( -0.378656, 0.178225, 2.41568));
    i_vertices.push_back( Vec3d( -0.391638, 0.263003, 2.07156));
    i_vertices.push_back( Vec3d( -0.405475, 0.327882, 1.72318));
    i_vertices.push_back( Vec3d( -0.413707, 0.40073, 1.21173));
    i_vertices.push_back( Vec3d( -0.382154, 0.474356, 0.701284));
    i_vertices.push_back( Vec3d( -0.30207, 0.569161, 0.19971));
    i_vertices.push_back( Vec3d( -0.203541, 0.684475, -0.294214));
    i_vertices.push_back( Vec3d( -0.105288, 0.81329, -0.784839));
    i_vertices.push_back( Vec3d( -0.0340096, 0.925021, -1.18432));
    i_vertices.push_back( Vec3d( 0.0264219, 1.04193, -1.58412));
    i_vertices.push_back( Vec3d( 0.0757033, 1.16426, -1.98382));
    i_vertices.push_back( Vec3d( 0.115412, 1.29266, -2.38268));
    i_vertices.push_back( Vec3d( 0.147909, 1.42739, -2.7801));
    i_vertices.push_back( Vec3d( 0.17159, 1.54614, -3.11535));
    i_vertices.push_back( Vec3d( 0.192911, 1.66858, -3.44944));
    i_vertices.push_back( Vec3d( 0.212798, 1.79392, -3.78255));
    i_vertices.push_back( Vec3d( 0.232007, 1.92103, -4.11501));


    std::cout << "num vertices: " << i_vertices.size() <<std::endl;
    // Enforce the total length
    Scalar length = 0.0;
    for ( int i = 0; i < nVertices - 1; i++ )
        length += ( i_vertices[i + 1] - i_vertices[i] ).norm();
    std::cout << "actual_length: " << length <<  std::endl;
    for ( int i = 0; i < nVertices; i++ )
        i_vertices[i] *= totalLength / length;
        
    VecXd dofs( nDOFs );
        
    // Initial strand position
    for ( int i = 0; i < dofs.size(); i += 4 )
        dofs.segment<3>( i ) = i_vertices[i / 4];
    
    DOFScriptingController* controller = new DOFScriptingController(  );
//    controller->freezeRootVertices<1>();
    controller->freezeVertices( nVertices - 1 );

    ElasticStrandParameters* params = new ElasticStrandParameters( 
                                            radiusA, 
                                            youngsModulus, 
                                            shearModulus, 
                                            density, 
                                            viscosity, 
                                            airDrag, 
                                            baseRotation );
        
    ElasticStrand* strand = new ElasticStrand( dofs, *params, controller );
    strand->setGlobalIndex( 0 );
    setRodCollisionParameters( *strand );
    Vec2Array kappas = strand->alterRestKappas();
    std::cout << "K_size " << kappas.size() <<std::endl;
    for( unsigned i = 0; i < kappas.size(); ++i ){
        strand->alterRestKappas()[i].setZero();
    }

    m_strands.push_back( strand );

    std::cout << "num strands = " << m_strands.size() <<'\n';
    std::cout << "num dofs per strand = " << nDOFs <<'\n';
    
    GravitationForce::setGravity( GetVecOpt("gravity").cast<Scalar>() );
}

void Knot::setupMeshes()
{}

bool Knot::executeScript()
{

     Vec3 zero(  0. , 0., 0. );
     Vec3 translate = GetVecOpt("translation");
     Mat3x id;
     id << 1,0,0,
           0,1,0,
           0,0,1;

    int vtx = 0; // GetIntOpt("nv") - 1;

      for(auto rd_itr = m_strands.begin(); rd_itr != m_strands.end(); ++ rd_itr)
      {
             SceneUtils::transformRodRootVtx( *rd_itr, id, zero, translate, vtx );
      }
    return true;
}

