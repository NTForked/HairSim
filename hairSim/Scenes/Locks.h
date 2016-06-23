#ifndef LOCKS
#define LOCKS

#include "../ProblemStepper.hh"

class Locks : public ProblemStepper
{
public:
    
    Locks();
    ~Locks();
    
protected:

	void layeredLockPacking( std::vector< std::vector< Vec3x > >& strands );
    void generateSinusoidalLocksVertices( std::vector< std::vector< Vec3x > >& strands );
	void includeHairTie();
    void setupStrands(); //TODO: virtual
    void setupMeshes(); //TODO: virtual
    bool executeScript(); // execute scripting and any per-step updating for problem
};
#endif 
