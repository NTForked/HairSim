#ifndef BRAID
#define BRAID

#include "../ProblemStepper.hh"

class Braid : public ProblemStepper
{
public:
    
    Braid();
    ~Braid();
    
protected:

    void generateBraidVertices( std::vector< std::vector< Vec3x > >& strands );
    void generateSinusoidalBraidVertices( std::vector< std::vector< Vec3x > >& strands );
    void loadNurbs();
	void includeHairTie();
    void setupStrands(); //TODO: virtual
    void setupMeshes(); //TODO: virtual
    bool executeScript(); // execute scripting and any per-step updating for problem
};
#endif 
