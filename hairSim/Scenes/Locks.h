#ifndef LOCKS
#define LOCKS

#include "Scene.h"

class Locks : public Scene
{
public:
    
    Locks();
    ~Locks();
    
protected:

	void layeredLockPacking( std::vector< std::vector< Vec3 > >& strands );
    void generateSinusoidalLocksVertices( std::vector< std::vector< Vec3 > >& strands );
	void includeHairTie();
    void setupStrands(); //TODO: virtual
    void setupMeshes(); //TODO: virtual
    bool executeScript(); // execute scripting and any per-step updating for problem
};
#endif 
