
#include "Render/Camera/ViewController.h"
#include "Render/Image.h"
#include "Scenes/Scene.h"
#include "Utils/Image.h"
#include "Utils/TypeDefs.h"

#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

// problems :
#include "Scenes/CousinIt.h"
#include "Scenes/PlayBack.h"
#include "Scenes/SingleContact.h"
#include "Scenes/Aleka.h"
#include "Scenes/Knot.h"
#include "Scenes/MultipleContact.h"
#include "Scenes/Braid.hh"
#include "Scenes/Twist.h"
#include "Scenes/Locks.h"
//

int g_problem_idx = -1;
int g_window_width = 512;
int g_window_height = 512;
int g_label = 0;
bool g_collisionThickness = true; // render with collision thickness or physical thickness

bool g_render = true;
bool g_paused = true;

bool g_restore_checkpoint = false;
bool g_dump_checkpoint = false;

time_t g_rawtime;
struct tm* g_timeinfo; // for time-stamped simulation capture & output directory

std::string g_outputdirectory; // directory to save to 
bool g_dumpcoord = false; // dump coords
bool g_dont_dumpmesh = true; // save space don't dump mesh

Scalar g_fps = 100; // framerate to generate movies &/or dumps with
int g_last_frame_num = -1; // last frame # that was output
int g_current_frame = 0;

Scene* g_ps;
ViewController controller;

using namespace std;

void createProblem()
{
    switch( g_problem_idx )
    {
        case 1:
            g_ps = new CousinIt();
            break;
        case 2:
            g_ps = new PlayBack();
            break;
        case 3:
            g_ps = new SingleContact();
            break;
        case 4:
            g_ps = new Aleka();
            break;
        case 5:
            g_ps = new Knot();
            break;
        case 6:
            g_ps = new MultipleContact();
            break;
        case 7:
            g_ps = new Braid();
            break;
        case 8:
            g_ps = new Twist();
            break;
        case 9:
            g_ps = new Locks();
            break;          
        default:
            cerr << "invalid Scene id" << endl;
            std::exit( EXIT_FAILURE );
            break;
    }
}

void setOptions()
{   // create options so they exist in scene files
    g_ps->AddOption( "render", "display output in OpenGL", g_render );
    g_ps->AddOption( "paused", "start simulation paused", g_paused );
    g_ps->AddOption( "dont_dump_mesh", "save space don't dump mesh data", g_dont_dumpmesh );
    g_ps->AddOption( "dump", "dump rod and mesh coordinates", g_dumpcoord );
    g_ps->AddOption( "fps", "frames per second for output movie and/or coordinate dump", g_fps );
}

void getOptions()
{   // read them back in case scene file changed them
    g_render = g_ps->GetBoolOpt("render");
    g_paused = g_render && g_ps->GetBoolOpt("paused");
    g_dont_dumpmesh = g_ps->GetBoolOpt("dont_dump_mesh");
    g_dumpcoord = g_ps->GetBoolOpt("dump");
    g_fps = g_ps->GetScalarOpt("fps");
}

void drawText( std::string s, int x, int y )
{
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0.0, g_window_width, 0.0, g_window_height);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glColor3f(0.9f, 0.9f, 0.9f);
    glRasterPos2i(x, g_window_height - y);
    void * font = GLUT_BITMAP_9_BY_15;
    for (std::string::iterator i = s.begin(); i != s.end(); ++i)
    {
        char c = *i;
        glutBitmapCharacter(font, c);
    }
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glEnable(GL_TEXTURE_2D);
}

void drawHUD()
{
    std::ostringstream oss; 
    if( nonLinearCallbackBogus) oss << "nonlinearCallback: ON";
    else oss << "nonlinearCallback: OFF";
    drawText( oss.str(), 7, 40 ); oss.str(""); oss.clear();

    if( hLoop ) oss << "hLoop: ON";
    else oss << "hLoop: OFF";
    drawText( oss.str(), 7, 55 ); oss.str(""); oss.clear();

    if( trackGeometricRelations ) oss << "trackGeometricRelations: ON";
    else oss << "trackGeometricRelations: OFF";
    drawText( oss.str(), 7, 70 ); oss.str(""); oss.clear();
}

void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glPushMatrix();
    
    controller.ApplyCamera();

    g_ps->render( g_window_width, g_window_height, g_label, g_collisionThickness );

    std::ostringstream oss; oss << "time: " << g_ps->getTime();
    drawText( oss.str(), 10, 20 );

    if( g_ps->isSimulated() ) drawHUD();

    glutSwapBuffers();
    glPopMatrix();
}

void storeCamera( void )
{
    std::ofstream out( ".viewer_state.txt" );

    Camera& c = controller.getCamera();

    Vec3d viewCenter = c.getViewCenter();
    out << viewCenter[0] << std::endl;
    out << viewCenter[1] << std::endl;
    out << viewCenter[2] << std::endl;

    Vec3d eye = c.getEye();
    out << eye[0] << std::endl;
    out << eye[1] << std::endl;
    out << eye[2] << std::endl;

    Vec3d up = c.getUp();
    out << up[0] << std::endl;
    out << up[1] << std::endl;
    out << up[2] << std::endl;

    int* w = &g_window_width;
    int* h = &g_window_height;
    c.getViewport( w, h );
    out << *w << std::endl;
    out << *h << std::endl;
}

void restoreCamera( void )
{
    std::ifstream in( ".viewer_state.txt" );

    if( !in.is_open() ) return;

    Camera& c = controller.getCamera();

    Vec3d viewCenter;
    in >> viewCenter[0];
    in >> viewCenter[1];
    in >> viewCenter[2];

    Vec3d eye;
    in >> eye[0];
    in >> eye[1];
    in >> eye[2];
    c.setEye( eye );

    Vec3d up;
    in >> up[0];
    in >> up[1];
    in >> up[2];
    c.setUp( up );

    c.setViewCenter( viewCenter );

    int w, h;
    in >> w;
    in >> h;

    g_window_width = w;
    g_window_height = h;
    
    controller.ApplyCamera();
    c.setPerspective(60, 1);
    const Scalar radius = controller.getBoundingRadius();
    c.setZClipping( 0.01 * radius, 3 * radius );
    c.setViewport( w, h );
    
    glutReshapeWindow( g_window_width, g_window_height );
    glutPostRedisplay();
}

void reshape( int w, int h )
{
    g_window_width = w;
    g_window_height = h;
    
    Camera& c = controller.getCamera();
    c.setPerspective(60, 1);
    const Scalar radius = controller.getBoundingRadius();
    c.setZClipping( 0.01 * radius, 3 * radius );
    c.setViewport( w, h );
    
    glutPostRedisplay();
}

void screenshot( void )
{
    static int index = 0;

    // get window width and height
    GLint view[4];
    glGetIntegerv( GL_VIEWPORT, view );
    int w = view[2];
    int h = view[3];

    // get pixels
    Image image( w, h );
    glReadPixels( 0, 0, w, h, GL_BGR, GL_FLOAT, &image(0,0) );

    stringstream filename;
    filename << "frames/frame" << setw(8) << setfill( '0' ) << index << ".tga";
    image.write( filename.str().c_str() );

    index++;
}

void output()
{
    if( ( g_dumpcoord || g_dump_checkpoint )
    {
        int steps_per_frame = -1;
        if( g_fps > 0 )
        {
            const double seconds_per_frame = 1.0 / ( (double) g_fps );
            const double steps_per_second = 1.0 / ( (double) g_ps->getDt() );
            steps_per_frame = floor( seconds_per_frame * steps_per_second + 0.5 ); // < rounds
        }
        
        const int frame = floor( g_ps->getTime() / g_ps->getDt() + 0.5 );
        
        if( frame % steps_per_frame == 0 && g_last_frame_num != frame )
        {
            g_last_frame_num = frame;
            mkdir( g_outputdirectory.c_str(), 0755 );
            const int file_width = 8;
            
            if ( g_dumpcoord )
            {
                g_ps->dumpRods( g_outputdirectory, g_current_frame, file_width );
                if (! g_dont_dumpmesh ){
                    g_ps->dumpMesh( g_outputdirectory, g_current_frame, file_width );
                }
                
                std::cout << "Saved coordinates of frame: " << g_current_frame << " @ time: " 
                    << g_ps->getTime() << " in directory " << g_outputdirectory << std::endl;
            }

            if( g_dump_checkpoint )
            {
                g_ps->checkpointSave( g_outputdirectory );
                g_ps->checkpointRestore( g_outputdirectory );
            }

            ++g_current_frame;
        }
    }
}

void stepOnce()
{
    g_ps->step();
    if( g_render )
    { // If in render mode, update the display
        glutPostRedisplay();
    }
    output();    
}

void idle()
{
    if( !g_paused )
    { // If the simulation isn't paused, take a timestep
        stepOnce();
    }
}

void setLighting()
{
    // Create a directional white light with a small ambient component
    glEnable( GL_LIGHT0 );
    GLfloat white_ambient[] = { 0.1, 0.1, 0.1, 1.0 };
    glLightfv( GL_LIGHT0, GL_AMBIENT, white_ambient );
    GLfloat white_diffuse[] = { 0.55, 0.55, 0.55, 1.0 };
    glLightfv( GL_LIGHT0, GL_DIFFUSE, white_diffuse );
    GLfloat upper_corner[] = { 1.0, 1.0, 1.0, 0.0 };
    glLightfv( GL_LIGHT0, GL_POSITION, upper_corner );
    
    // Create a much weaker direction light
    glEnable( GL_LIGHT1 );
    GLfloat weak_white_diffuse[] = { 0.3, 0.3, 0.3, 1.0 };
    glLightfv( GL_LIGHT1, GL_DIFFUSE, weak_white_diffuse );
    GLfloat negative_z[] = { 0.0, 0.0, 1.0, 0.0 };
    glLightfv( GL_LIGHT1, GL_POSITION, negative_z );
    
    glShadeModel( GL_FLAT );
}

void centerObject()
{
    Vec3d center = Vec3d::Zero();
    g_ps->getCenter( center );
    controller.setCenterMode( ViewController::CENTER_OBJECT );
    controller.setViewCenter( center );
    
    Scalar radius = 0.0;
    g_ps->getRadius( radius, center );
    if( radius == 0.0 ) radius = 0.1;
    controller.setBoundingRadius( radius );
}

void initCamera()
{
    controller.setViewDirection( Vec3d( 0, 0, -2 ) );
    centerObject();
}

void scaleMousePos( int x, int y, Scalar& xx, Scalar& yy )
{
    int w, h;
    controller.getCamera().getViewport( &w, &h );
    
    xx = 2 * x / (Scalar) (w - 1) - 1.0;
    yy = 2 * (h - y - 1) / (Scalar) (h - 1) - 1.0;
}

void mouse( int button, int state, int x, int y )
{
    const bool zooming = (button == GLUT_MIDDLE_BUTTON) || ((button == GLUT_LEFT_BUTTON) && (glutGetModifiers() & GLUT_ACTIVE_CTRL));
    const bool translating = (button == GLUT_LEFT_BUTTON) && (glutGetModifiers() & GLUT_ACTIVE_SHIFT);
    const bool scripting = (button == GLUT_LEFT_BUTTON) && (glutGetModifiers() & GLUT_ACTIVE_ALT);
    const bool rotating = (button == GLUT_LEFT_BUTTON) && (glutGetModifiers() == 0);
    
    scalar xx, yy;
    scaleMousePos(x, y, xx, yy);
    if( state == GLUT_DOWN ){
        if( translating ) controller.beginTranslationDrag(xx, yy);
        if( zooming ) controller.beginZoomDrag(xx, yy);
        if( rotating ) controller.beginRotationDrag(xx, yy);
    }
    else{
        controller.endTranslationDrag( xx, yy );
        controller.endRotationDrag( xx, yy );
        controller.endZoomDrag( xx, yy );
    }
    
    glutPostRedisplay();
}

void motion( int x, int y )
{
    Scalar xx, yy;
    scaleMousePos( x, y, xx, yy );
    controller.updateDrag( xx, yy );
    glutPostRedisplay();
}

void SpecialInput(int key, int x, int y)
{
    switch(key)
    {
        case GLUT_KEY_LEFT:
        case GLUT_KEY_RIGHT:
        case GLUT_KEY_DOWN:
        case GLUT_KEY_UP:
            break;
    }
}

void keyboard( unsigned char key, int, int )
{
    switch( key )
    {
        case 'q':
        case 27: // ESC key
            std:exit(EXIT_SUCCESS);
            break;
        case 'p':
        case ' ':
            g_paused = !g_paused;
            break;
        case 't':
            g_collisionThickness = !g_collisionThickness;
            break;
        case 'f':
            stepOnce();
            break;
        case 'm':
            screenshot();
            break;            
        case 'u':
            if( g_render )
            {
                g_label += 1;
                glutPostRedisplay();
            }
            break; 
        case 'w':
            storeCamera();
            break;  
        case 'e':
            restoreCamera();
            break;
    }
}

void initializeOpenGL( int argc, char** argv )
{
    glutInit( &argc, argv );
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_ALPHA);
    glutInitDisplayString("rgba double depth alpha samples>=16");
    glutInitWindowPosition( 50, 50);
    glutInitWindowSize( g_window_width, g_window_height );
    glutCreateWindow( argv[0] );
    // glutFullScreen();
    glEnable( GL_DEPTH_TEST );
    glDepthFunc( GL_LESS );
    glClearColor(161.0/255.0, 161.0/255.0, 161.0/255.0, 0.0);
    // glClearColor(0.0/255.0,0.0/255.0,0.0/255.0,0.0);
    // glClearColor( 0.15, 0.15, 0.15, 0.0 );
    
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutKeyboardFunc(keyboard);
    glutIdleFunc(idle);
    glutSpecialFunc(SpecialInput);
}

void printStamp()
{
    time_t g_rawEndtime;
    time( &g_rawEndtime );
    std::cout << "# wallClock seconds passed: " << difftime( g_rawEndtime, g_rawtime ) << std::endl;
}

void cleanup()
{
    printStamp();
    delete g_ps;
}

void createOptionsFile( const std::string& filename )
{
    std::ofstream file;
    file.open( filename.c_str() );
    if (!file.is_open())
    {
        std::cerr << "Failed to open file " << filename << std::endl;
        return;
    }
    
    g_ps->PrintOptions(file);
    
    file.close();
    std::cout << "Generated options file for " << g_ps->m_problemName << ": " << filename << std::endl;
}

int parseCommandLine( int argc, char** argv )
{    
    try
    {
        TCLAP::CmdLine cmd("hairSim");
        TCLAP::ValueArg<int> run( "r", "run", "Run a problem", true, -1, cmd );
        TCLAP::ValueArg<std::string> file( "f", "file", "Options file for a problem", true, "", "string", cmd );
        TCLAP::ValueArg<bool> dumpcoord ( "d", "dumpcoord", "Dump coordinates of all rods and meshes at each frame", false, false, "boolean", cmd );
        cmd.parse(argc, argv);
        
        if( run.isSet() )
        {
            g_problem_idx = run.getValue();
            createProblem();
            setOptions();            
            createOptionsFile( file.getValue() );
            if( g_ps->LoadOptions(file.getValue()) == -1 )
            {
                return -1;
            }
            getOptions();
            
            if( dumpcoord.isSet() )
            {
                g_dumpcoord = true;
                
                // dump copy of options used for the record
                mkdir( g_outputdirectory.c_str(), 0755 );
                std::stringstream name;
                name << g_outputdirectory << "/" << "config_file";
                std::ofstream file;
                file.open(name.str().c_str());
                if (!file.is_open())
                {
                    std::cerr << "Failed to open file " << name.str() << std::endl;
                    return -1;
                }
                
                file << "# PROBLEM : " << g_ps->m_problemName << std::endl << "# SAVED TO : " << g_outputdirectory.c_str() <<  std::endl;
                createOptionsFile( file );
                file.close();
                
                std::cout << "Dumped options for problem" << g_ps->m_problemName << " to file " << name.str() << std::endl;
            }
            
            return idx;
        }

        std::cerr << cmd.getProgramName() << ": missing operand" << std::endl 
        << "Try `" << cmd.getProgramName() << " --help'" << " for more information" << std::endl;
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "ERROR: " << e.argId() << std::endl << "       " << e.error() << std::endl;
    }
    
    return -1;
}

void printTimeStamp()
{
    std::cout << std::setfill('0') << (1900 + g_timeinfo->tm_year) << "/" << std::setw(2) << (1 + g_timeinfo->tm_mon) << "/";
    std::cout << (g_timeinfo->tm_mday) << ", " << std::setw(2) << (g_timeinfo->tm_hour) << ":" << std::setw(2) << (g_timeinfo->tm_min);
    std::cout << ":" << std::setw(2) << (g_timeinfo->tm_sec) << std::endl;    
}

void printCommandLineSplashScreen()
{
    #ifdef DEBUG
        std::cout << "# build mode: DEBUG" << std::endl;
    #else
        std::cout << "# build mode: RELEASE" << std::endl;
    #endif
    std::cout << "# timestamp: "; printTimeStamp();
    std::cout << "# problem: " << g_ps->m_problemName << std::endl;
}

std::string generateOutputDirName()
{ 
    time( &g_rawtime );
    g_timeinfo = localtime( &g_rawtime );
    assert( g_timeinfo != NULL );
    std::stringstream datestream;
    datestream.fill('0');
    datestream << std::setw(4) << (1900 + g_timeinfo->tm_year) << "_" << std::setw(2) << (1 + g_timeinfo->tm_mon) << "_";
    datestream << std::setw(2) << (g_timeinfo->tm_mday) << "_" << std::setw(2) << (g_timeinfo->tm_hour) << "_" << std::setw(2) << (g_timeinfo->tm_min) << "_";
    datestream << std::setw(2) << (g_timeinfo->tm_sec) << "_" << "output";
    return datestream.str();
}

int main( int argc, char** argv )
{
    // make time-stamped directory name
    g_outputdirectory = generateOutputDirName();
    
    atexit( cleanup );
    
    if( parseCommandLine( argc, argv ) < 0 )
    {
        return -1;
    }
    
    g_ps->setup();
    if( g_restore_checkpoint ) g_ps->checkpointRestore();
    
    printCommandLineSplashScreen();
    
    if( g_render )
    {
        initializeOpenGL( argc, argv );
        initCamera();
        setLighting();
        glutMainLoop();
    }
    else{
        g_paused = false;
        for(;;) idle();
    }
 
    return 0;
}
