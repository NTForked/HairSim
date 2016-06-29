#ifndef MULTIPLE_CONTACT
#define MULTIPLE_CONTACT

#include "Scene.h"

class MultipleContact : public Scene
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
