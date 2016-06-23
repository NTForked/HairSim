#ifndef SINGLE_CONTACT
#define SINGLE_CONTACT

#include "../ProblemStepper.hh"

class SingleContact : public ProblemStepper
{
// One strand falling orthogonally on
// another strand fixed at both ends
public:
    SingleContact();
    ~SingleContact();
    
protected:
    Scalar m_radius;
    void setupStrands();
    void setupMeshes();
    bool executeScript();
};

#endif
