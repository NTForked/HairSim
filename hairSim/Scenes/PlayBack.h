#ifndef PLAYBACK
#define PLAYBACK

#include "Scene.h"

class Playback : public Scene
{
public:
    
    Playback();
    ~Playback();
    
protected:

    void setupStrands(); 
    void setupMeshes(); 
    bool executeScript(); // execute scripting and any per-step updating for problem

private:
    int m_current_frame;
    void loadRods(int frame);
    void loadMeshes(int frame);

};
#endif 
