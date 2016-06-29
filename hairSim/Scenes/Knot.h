#ifndef KNOT
#define KNOT

#include "Scene.h"

class Knot : public Scene
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
