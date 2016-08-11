#include <iostream>
#include <math.h>
#include <string>
#include <GLUT/GLUT.h>
#include "bmphandler.hpp"
using namespace std;
GLint universe_time = 0;
string imagePath = "/Users/anakin/github/network-summer/BC-summer/openGL/openGL/pic";
GLuint textureIndex[11];
GLUquadric *quad;
string matterNames[11]  = {"sun.bmp", "mercury.bmp", "venus.bmp", "earth.bmp", "mars.bmp", "jupiter.bmp", "saturn.bmp", "uranus.bmp", "neptune.bmp", "moon.bmp", "star.bmp"};
GLfloat period[10] = {0.0f,100.0f,200.0f,300.0f,350.0f,500.0f,650.0f,900.0f,1000.0f,300.0f};
GLfloat radius[10] = {0.0f,2.6f,3.6f,5.1f,6.6f,8.6f,10.6f,12.6f,14.6f,1.0f};
GLfloat sphereRadius[10] = {0.0f,0.2f,0.4f,0.5f,0.35f,0.65f,0.75f,0.3f,0.35f,0.2f};
float fTranslate;
float fRotate;
float fScale     = 1.0f;	// set inital scale value to 1.0f

bool bPersp = false;
bool bAnim = false;
bool bWire = false;

int wHeight = 0;
int wWidth = 0;
float teapotX = 0;
float teapotY = 0;
float moveSize = 0.5;
float teapotAng = 0;

//todo
//hint: some additional parameters may needed here when you operate the teapot
GLfloat calAngle(GLfloat thisperiod)
{
   GLfloat res;
   res = universe_time / thisperiod;
   res = (res-(int)res)*360;
   return res;
}

void Draw_Scene()
{
   //sun
   glPushMatrix();
   glRotatef(calAngle(200), 0.0f, 0.0f, 0.0f);
   glBindTexture(GL_TEXTURE_2D, textureIndex[0]);
   gluSphere(quad, 1.5, 48, 24);
   glPopMatrix();
   
   
   //planets
   for(int i = 1; i < 9; i++){
       glPushMatrix();
       GLfloat angle = calAngle(period[i]);
       glRotatef(angle, 0.0f, 0.0f, 1.0f);
       glTranslatef(radius[i], 0, 0);
       glRotatef(calAngle(200), 0.0f, 0.0f, 0.0f);
       glBindTexture(GL_TEXTURE_2D, textureIndex[i]);
       gluSphere(quad, sphereRadius[i], 48, 24);
       glPopMatrix();
   }
   //moon
   glPushMatrix();
   glRotatef(calAngle(period[3]), 0.0f, 0.0f, 1.0f);
   glTranslatef(radius[3], 0, 0);
   GLfloat angle = calAngle(period[9]);
   glRotatef(angle, 0.0f, 0.0f, 1.0f);
   glTranslatef(radius[9], 0, 0);
   glRotatef(calAngle(200), 0.0f, 0.0f, 0.0f);
   glBindTexture(GL_TEXTURE_2D, textureIndex[9]);
   gluSphere(quad, sphereRadius[9], 48, 24);
   glPopMatrix();
   
   //ellipse
   glPushMatrix();
   angle = calAngle(100);
   glRotatef(angle, 0.0f, 0.0f, 1.0f);
   GLfloat xx = 15*cos(2*angle*3.1415926/720);
   GLfloat yy = 9*sin(2*angle*3.1415926/720);
   glTranslatef(sqrtf(xx*xx+yy*yy), 0, 0);
   glBindTexture(GL_TEXTURE_2D, textureIndex[10]);
   gluSphere(quad, 0.1, 48, 24);    
   glPopMatrix();
   
   
   
   //radius
   for (int j = 1; j < 9; j++) {
       glBegin(GL_LINE_LOOP);
       //glVertex2f(0, 0);
       float x[40],y[40];
       for(int i = 1; i < 40; i++)
       {
           x[i] = radius[j]*cos(2*i*3.1415926/40);
           y[i] = radius[j]*sin(2*i*3.1415926/40);
           glVertex2f(x[i], y[i]);
       }
       glEnd();

   }
   glPushMatrix();
   glRotatef(calAngle(period[3]), 0.0f, 0.0f, 1.0f);
   glTranslatef(radius[3], 0, 0);
   glBegin(GL_LINE_LOOP);
   //glVertex2f(0, 0);
   float x[40],y[40];
   for(int i = 1; i < 40; i++)
   {
       x[i] = radius[9]*cos(2*i*3.1415926/40);
       y[i] = radius[9]*sin(2*i*3.1415926/40);
       glVertex2f(x[i], y[i]);
   }
   glEnd();
   glPopMatrix();

}

void updateView(int width, int height)
{
   glViewport(0,0,width,height);						// Reset The Current Viewport

   glMatrixMode(GL_PROJECTION);						// Select The Projection Matrix
   glLoadIdentity();									// Reset The Projection Matrix

   float whRatio = (GLfloat)width/(GLfloat)height;

   if (bPersp){
       gluPerspective(30, whRatio, 1, 100);
       //todo when 'p' operation, hint: use FUNCTION gluPerspective
   }
   else
       glOrtho(-3 ,3, -3, 3,-100,100);

   glMatrixMode(GL_MODELVIEW);							// Select The Modelview Matrix
}

void reshape(int width, int height)
{
   if (height==0)										// Prevent A Divide By Zero By
   {
       height=1;										// Making Height Equal One
   }

   wHeight = height;
   wWidth = width;

   updateView(wHeight, wWidth);
}

//void idle()
//{
//    glutPostRedisplay();
//}

float eye[] = {0, 0, 8};
float center[] = {0, 0, 0};
//todo; hint: you may need another ARRAY when you operate the teapot

void key(unsigned char k, int x, int y)
{
   switch(k)
   {
       case 27:
       case 'q': {exit(0); break; }
       case 'p': {bPersp = !bPersp; updateView(wHeight, wWidth);break; }

       case ' ': {bAnim = !bAnim; break;}
       //case 'o': {bWire = !bWire; break;}

       //the direction
       case 'a': {//todo, hint: eye[] and center[] are the keys to solve the problems
           eye[0]+=moveSize;
           //center[0]+=moveSize;
           break;
       }
       case 'd': {//todo
           eye[0]-=moveSize;
           //center[0]-=moveSize;
           break;
       }
       case 'w': {//todo
           eye[1]-=moveSize;
           //center[1] -= moveSize;

           break;
       }
       case 's': {//todo
           eye[1]+=moveSize;
           //center[1]+=moveSize;

           break;
       }


       case 'z': {//todo
           eye[2]+=moveSize;
           center[2] += moveSize;
           break;
       }
       case 'c': {//todo
           eye[2]-=moveSize;
           center[2] -= moveSize;
           break;
       }

           //deal with the teapot
       case 'l': {//todo, hint:use the ARRAY that you defined, and notice the teapot can NOT be moved out the range of the table.
           if(teapotX <= 2.5)
               teapotX+=moveSize;
           break;
       }
       case 'j': {//todo
           if(teapotX >= -2.5)
               teapotX-=moveSize;
           break;
       }
       case 'i': {//todo
           if(teapotY <= 2.5)
               teapotY+=moveSize;
           break;
       }
       case 'k': {//todo
           if(teapotY >= -2.5)
               teapotY-=moveSize;
           break;
       }
       case 'e': {//todo
           teapotAng+=10;
           break;
       }
   }
}


void redraw()
{

   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   glEnable(GL_TEXTURE_2D);
   glEnable(GL_DEPTH_TEST);
   glEnable(GL_LIGHTING);
   glLoadIdentity();									// Reset The Current Modelview Matrix

   gluLookAt(eye[0], eye[1], eye[2],
             center[0], center[1], center[2],
             0, 1, 0);

//    if (bWire) {
//        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
//    }
//    else {
//        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
//    }

   
   GLfloat white[] = { 1.0, 1.0, 1.0, 1.0 };
   GLfloat red[] = { 1.0, 0.8, 0.8, 1.0 };
   //GLfloat light_pos[] = {5,5,5,1};
   GLfloat light_pos[] = {0,0,0,1};
   GLfloat light_pos2[] = {5,5,5,1};

   glLightfv(GL_LIGHT0, GL_POSITION, light_pos);
   glLightfv(GL_LIGHT0, GL_AMBIENT, red);
   glEnable(GL_LIGHT0);

   glLightfv(GL_LIGHT1, GL_POSITION, light_pos2);
   glLightfv(GL_LIGHT1, GL_AMBIENT, white);
   glEnable(GL_LIGHT1);

   //	glTranslatef(0.0f, 0.0f,-6.0f);			// Place the triangle at Center
   glRotatef(fRotate, 0, 1.0f, 0);			// Rotate around Y axis
   //glRotatef(-90, 1, 0, 0);
   glScalef(0.2, 0.2, 0.2);
   Draw_Scene();						// Draw Scene
   
   

   if (bAnim) fRotate    += 0.5f;

   //todo; hint: when you want to rotate the teapot, you may like to add another line here =)
   glutSwapBuffers();
}
void setTexture()
{
   glGenTextures(10, textureIndex);
   for(int i = 0; i < 11; i++){
       string bmpPic;
       RGBpixmap bmpPicHandler;
       bmpPic = imagePath + "/" + matterNames[i];
       bmpPicHandler.readBMPFile(bmpPic);
       bmpPicHandler.setTexture(textureIndex[i]);
   }
   
   quad = gluNewQuadric();
   gluQuadricTexture(quad, GL_TRUE);
   glEnable(GL_TEXTURE_2D);
}
void initialize()
{
   //set the background color
   glEnable(GL_DEPTH_TEST);
   glClearColor(0, 0, 0, 0);
   glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
   setTexture();
}
void Timer(int para)
{
   universe_time += 1;
   glutPostRedisplay();
   glutTimerFunc(40, Timer, 1);
}

int main (int argc,  char *argv[])
{
   glutInit(&argc, argv);
   glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
   glutInitWindowSize(800,800);
   int windowHandle = glutCreateWindow("Simple GLUT App");
   initialize();
   
   glutDisplayFunc(redraw);
   glutReshapeFunc(reshape);
   glutKeyboardFunc(key);
   //glutIdleFunc(idle);
   glutTimerFunc(40, Timer, 1);

   glutMainLoop();
   return 0;
}


