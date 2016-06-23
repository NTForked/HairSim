#ifndef KNOT
#define KNOT

#include "../ProblemStepper.hh"

class Knot : public ProblemStepper
{
// One strand knotted and pulled tight
public:
    Knot();
    ~Knot();
    
protected:
    Scalar m_radius;
    void setupStrands();
    void setupMeshes();
    bool executeScript();
};

#endif
