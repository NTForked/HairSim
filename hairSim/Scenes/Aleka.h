#ifndef ALEKA
#define ALEKA

#include "Scene.h"

class Aleka : public Scene
{
// M strands passing over
// N strands fixed at both ends
public:
    Aleka();
    ~Aleka();
    
protected:
    Scalar m_radius;
    void setupStrands();
    void setupMeshes();
    bool executeScript();
};

#endif
