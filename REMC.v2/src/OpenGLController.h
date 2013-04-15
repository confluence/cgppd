// OpenGLController.h
// this file contains all global functions for displaying whats going on in opengl,
// its non essential for the contant of this code so its all been moved here to keep
// it separate from the other stuff

#ifndef OPENGL_H
#define OPENGL_H

#include <GL/glut.h>
#include <GL/gl.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include "Camera.h"
#include <vector>
#include <ctime>
#include <linux/kernel.h>
//#include <linux/time.h>

#define WIDTH           800
#define HEIGHT          800
#define MOVESPEED  0.05f  //movement
#define MSPEED     0.02f  //mouse
#define MACC       0.8f  //mouse acceleration

#define M_ROTATE        1001
#define M_TRANSLATE     1002
#define M_SCALE         1003
#define M_QUIT          1004
#define M_COLOR         1005
#define M_ADD           1006
#define M_REMOVE        1007
#define M_SAVE          1008
#define M_LOAD          1009
#define M_DESELECT      1010
#define M_SELECT        1011
#define M_WIREFRAME     1012
#define M_FILLED        1013
#define M_NORMALS       1014
#define M_AXES          1015
#define M_INVSELECTION  1016

#define AXISDIM 1500

// TODO replace with pointer attribute
// extern Replica replica[REPLICA_COUNT];

class OpenGLController
{
public:
    bool displayminimum;

    Camera camera;

    Replica (*replica)[REPLICA_COUNT]; // TODO is this right?
    Replica * GLreplica;

    //modes
    bool b_mlook;
    bool b_rotate;
    bool b_translate;
    bool b_scale;

    //draw options
    bool drawAxes;
    bool drawNormals;
    bool drawFilled;
    bool drawWire;
    bool pulse;
    bool glow;
    bool fog;

    int mousex;
    int mousey;
    int windowh;
    int windoww;
    bool showSelection;
    int gl_replicaCount;
    int gl_boundingValue;
    int dr; //replica to display
    // input control variables 1,0,-1 indicate states
    int moveFrontback;
    int moveLeftRight;
    float deltaRotate;
    unsigned char lastkey;
    int glutWindowID;

    float LightPos0[];			// Light Position
    float LightPos1[];			// Light Position
    float LightPos2[];			// Light Position
    float materialShininess;										// Material shininess

    float COLOURS[5][3];

    char *windowtitle;
    // initialization function

    OpenGLController();
    ~OpenGLController();

    void init(int argc, char **argv, argdata parameters);

    void GlutInit(int windowheight,int windowwidth, char* title);
    void GlutDisplay();
    void GlutResize(int width, int height);
    void GlutMouseFunction(int btn, int state, int x, int y);
    void GlutIdleFunction();
    void GLUTMouseActiveMotion(int x, int y);
    void GLUTMousePassiveMotion(int x, int y);
    void GLUTMouse(int btn, int state, int x, int y);
    void GLUTKeyboardPress(unsigned char x , int y ,int z);
    void GLUTKeyboardUp(unsigned char x , int y ,int z);
    void GLUTKeyboardSpecial(int x , int y ,int z);
    void GLUTcreateMenu();
    void GLUTmenu(int);
    void GLUTdrawaxes();
}

#endif
