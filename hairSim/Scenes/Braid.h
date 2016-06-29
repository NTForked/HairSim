#ifndef BRAID
#define BRAID

#include "Scene.h"

class Braid : public Scene
{
public:
    
    Braid();
    ~Braid();
    
protected:

    void generateBraidVertices( std::vector< std::vector< Vec3 > >& strands );
    void generateSinusoidalBraidVertices( std::vector< std::vector< Vec3 > >& strands );
    void loadNurbs();
	void includeHairTie();
    void setupStrands(); //TODO: virtual
    void setupMeshes(); //TODO: virtual
    bool executeScript(); // execute scripting and any per-step updating for problem
};
#endif 
