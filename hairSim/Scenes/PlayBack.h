#ifndef PLAYBACK
#define PLAYBACK

#include "Scene.h"

class PlayBack : public Scene
{
public:
    
    PlayBack();
    ~PlayBack();
    
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
