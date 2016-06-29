#ifndef SINGLE_CONTACT
#define SINGLE_CONTACT

#include "Scene.h"

class SingleContact : public Scene
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
