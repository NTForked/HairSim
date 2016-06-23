#ifndef MULTIPLE_CONTACT
#define MULTIPLE_CONTACT

#include "../ProblemStepper.hh"

class MultipleContact : public ProblemStepper
{
// One strand falling orthogonally on
// another strand fixed at both ends
public:
    MultipleContact();
    ~MultipleContact();
    
protected:
    Scalar m_radius;
    void setupStrands();
    void setupMeshes();
    bool executeScript();
};

#endif
