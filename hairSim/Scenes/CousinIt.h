#ifndef COUSINIT
#define COUSINIT

#include "../ProblemStepper.hh"

class CousinIt : public ProblemStepper
{
public:
    
    CousinIt();
    ~CousinIt();
    
protected:

    void generateNormalSamples( Scalar hair_region, int num_hairs, std::vector<Vec3x>& normals );
    void loadStrands();
    void setupStrands(); //TODO: virtual
    void setupMeshes(); //TODO: virtual
    bool executeScript(); // execute scripting and any per-step updating for problem
};
#endif 
