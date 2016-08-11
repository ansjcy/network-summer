//// glutEx1.cpp : ¶¨Òå¿ØÖÆÌ¨Ó¦ÓÃ³ÌÐòµÄÈë¿Úµã¡£
////
//#include <stdlib.h>
//#include <stdio.h>
//#include <string.h>
//#include <math.h>
//#include <GLUT/GLUT.h>
//#include <OpenGL/OpenGL.h>
//
//int wHeight = 0;
//int wWidth = 0;
//
//void startPicking(int cursorX, int cursorY);
//void stopPicking();
//bool mypickup[3][3] = { false, false, false, false, false, false, false, false, false };
//bool mypickdown[3][3] = { false, false, false, false, false, false, false, false, false };
//float wmat_diffuse[] = { 0.0, 0.0, 0.0 };
//float wmat_specular[] = { 2.0, 1.0, 1.0, 1.0 };
//
//
//void glCircle3i(GLint x, GLint y, GLint radius) {
//    float angle;
//    glLineWidth(1.0f);
//    glBegin(GL_POLYGON);
//    int i = 0;
//    for(i = 0; i < 100; i++) {
//        angle = i*2*M_PI/100;
//        glVertex2f(x + (cos(angle) * radius), y + (sin(angle) * radius));
//    }
//    glEnd();
//}
//
//void drawCircles()
//{
////    for (i = 0; i < size; i++){
////        glPushMatrix();
////        int posX = n[i].position.x;
////        int posY = n[i].position.y;
////        
////        glColor3f(1.0f, 0.2f, 0.2f);
////        glCircle3i(posX, posY, 10);
////        glColor3f(0.1f, 0.1f, 0.5f);
////        glPopMatrix();
////    }
//    
//    glPushMatrix();
//
//    glColor3f(1.0f, 0.2f, 0.2f);
//    glRotatef(90, 1, 0, 0);
//    glCircle3i(0, 0, 10);
//    glColor3f(0.1f, 0.1f, 0.5f);
//    glPopMatrix();
//    
//}
//
////void Draw_Desk();
//
//void drawThings() // This function draws a triangle with RGB colors
//{
//    for (int i = 0; i<1; i++)
//    for (int j = 0; j<1; j++)
//    {
//        
//        glPushName(i * 3 + j);
//        if (mypickdown[i][j] == true)
//        {
////            glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
////            glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
//            //glMateriali(GL_FRONT, GL_SHININESS, 20);
//            glPushMatrix();
//            glTranslatef(-1 + i, -1 + j, 5.5);
//            glRotatef(90, 1, 0, 0);
//            glColor3f(1.0, 1.0, 0.0);
//            //drawCircles();
//            glCircle3i(0, 0, 10);
//            glPopMatrix();
//        }
//        else
//        {
////            glMaterialfv(GL_FRONT, GL_DIFFUSE, wmat_diffuse);
////            glMaterialfv(GL_FRONT, GL_SPECULAR, wmat_specular);
//            //glMateriali(GL_FRONT, GL_SHININESS, 20);
//            glPushMatrix();
//            glTranslatef(-1 + i, -1 + j, 5.5);
//            glRotatef(90, 1, 0, 0);
//            glColor3f(1.0, 1.0, 1.0);
//            //drawCircles();
//            glCircle3i(0, 0, 10);
//            glPopMatrix();
//        }
//        glPopName();
//    }
////    for (int i = 0; i<3; i++)
////    for (int j = 0; j<3; j++)
////    {
////        glPushName(9 + i * 3 + j);
////        if (mypickup[i][j] == true)
////        {
////            glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
////            glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
////            glMateriali(GL_FRONT, GL_SHININESS, 70);
////            glPushMatrix();glColor3f(1.0, 0, 0);
////            glTranslatef(-1 + i, -1 + j, 7.5);
////            drawCircles();
////            glPopMatrix();
////        }
////        else
////        {
////            glMaterialfv(GL_FRONT, GL_DIFFUSE, wmat_diffuse);
////            glMaterialfv(GL_FRONT, GL_SPECULAR, wmat_specular);
////            glMateriali(GL_FRONT, GL_SHININESS, 70);
////            glPushMatrix();
////            glColor3f(1.0, 1.0, 1.0);
////            glTranslatef(-1 + i, -1 + j, 7.5);
////            drawCircles();
////            glPopMatrix();
////        }
////        glPopName();
////        
////    }
//
//}
//
//void updateView(int width, int height)
//{
//    glViewport(0,0,width,height);						// Reset The Current Viewport
//    
//    glMatrixMode(GL_PROJECTION);						// Select The Projection Matrix
//    glLoadIdentity();									// Reset The Projection Matrix
//    
////    float whRatio = (GLfloat)width/(GLfloat)height;
////    if (bPersp) {
////        gluPerspective(45.0f, whRatio,0.1f,100.0f);
////        //glFrustum(-3, 3, -3, 3, 3,100);
////    } else {
//        glOrtho(-3 ,3, -3, 3,-100,100);
////    }
//    
//    glMatrixMode(GL_MODELVIEW);							// Select The Modelview Matrix
//}
//
//void reshape(int width, int height)
//{
//    if (height==0)										// Prevent A Divide By Zero By
//    {
//        height=1;										// Making Height Equal One
//    }
//    
//    wHeight = height;
//    wWidth = width;
//    
//    updateView(wHeight, wWidth);
//}
//
//void idle()
//{
//    glutPostRedisplay();
//}
//
//float center[] = {0, 0, 0};
//float eye[] = {0, 0, 20};
//
//void key(unsigned char k, int x, int y)
//{
//    switch(k)
//    {
//            case 27:
//            case 'q': {exit(0); break; }
//            case '0': {}
//            
//            case 'a': {
//                eye[0] -= 0.2f;
//                center[0] -= 0.2f;
//                break;
//            }
//            case 'd': {
//                eye[0] += 0.2f;
//                center[0] += 0.2f;
//                break;
//            }
//            case 'w': {
//                eye[1] -= 0.2f;
//                center[1] -= 0.2f;
//                break;
//            }
//            case 's': {
//                eye[1] += 0.2f;
//                center[1] += 0.2f;
//                break;
//            }
//            case 'z': {
//                eye[2] -= 0.2f;
//                center[2] -= 0.2f;
//                break;
//            }
//            case 'c': {
//                eye[2] += 0.2f;
//                center[2] += 0.2f;
//                break;
//            }
//    }
//    
//    updateView(wHeight, wWidth);
//}
//
//void redraw()
//{
//    
//    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//    glLoadIdentity();									// Reset The Current Modelview Matrix
//    
//    gluLookAt(eye[0], eye[1], eye[2],
//              center[0], center[1], center[2],
//              0, 1, 0);				// ³¡¾°£¨0£¬0£¬0£©µÄÊÓµãÖÐÐÄ (0,5,50)£¬YÖáÏòÉÏ
//    
//    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
//    
//    //	glTranslatef(0.0f, 0.0f,-6.0f);			// Place the triangle at Center
//    //glRotatef(fRotate, 0, 1.0f, 0);			// Rotate around Y axis
//    glRotatef(-90, 1, 0, 0);
//    glScalef(0.2, 0.2, 0.2);
//    drawThings();						// Draw triangle
//    
//    glutSwapBuffers();
//}
//#define BUFSIZE 512
//GLuint selectBuf[BUFSIZE];
//
//void startPicking(int cursorX, int cursorY)
//{
//    printf("cx=%d, cy=%d\n", cursorX, cursorY);
//    GLint viewport[4];
//    glSelectBuffer(BUFSIZE, selectBuf);
//    glRenderMode(GL_SELECT);
//    glMatrixMode(GL_PROJECTION);
//    glPushMatrix();
//    glLoadIdentity();
//    glGetIntegerv(GL_VIEWPORT, viewport);
//    printf("viewport:: %d, %d, %d, %d\n", viewport[0], viewport[1], viewport[2], viewport[3]);
//    gluPickMatrix(cursorX, viewport[3] - cursorY, 1, 1, viewport);
//    gluPerspective(60, 1, 0.1, 1000);
//    glMatrixMode(GL_MODELVIEW);
//    glInitNames();
//    drawThings();
//    stopPicking();
//}
//void processHits(GLint hits, GLuint buffer[])
//{
//    GLuint names, *ptr, minZ, *ptrNames, numberOfNames;
//    
//    printf("hits = %d\n", hits);
//    ptr = (GLuint *)buffer;
//    minZ = 0xffffffff;
//    int i, j, m;
//    
//    for (i = 0; i<hits; i++) {
//        names = *ptr;
//        ptr++;
//        
//        if (*ptr < minZ) {
//            numberOfNames = names;
//            minZ = *ptr;
//            ptrNames = ptr + 2;
//        }
//        
//        ptr += names + 2;
//    }
//    
//    printf("The closest are ");
//    ptr = ptrNames;
//    for (j = 0; j<numberOfNames; j++, ptr++) {
//        printf("%d ", *ptr);
//        m = *ptr;
//    }
//    printf("\n");
//    if (m<9)
//    {
//        i = m / 3;
//        j = m % 3;
//        mypickdown[i][j] = !mypickdown[i][j];
//    }
//    else if (m >= 9 && m <= 17)
//    {
//        i = (m - 9) / 3;
//        j = (m - 9) % 3;
//        
//        mypickup[i][j] = !mypickup[i][j];
//    }
//}
//
//void stopPicking() {
//    int hits;
//    
//    // restoring the original projection matrix
//    glMatrixMode(GL_PROJECTION);
//    glPopMatrix();
//    glMatrixMode(GL_MODELVIEW);
//    glFlush();
//    
//    // returning to normal rendering model
//    hits = glRenderMode(GL_RENDER);
//    
//    // if there are hits process them
//    if (hits != 0)
//    processHits(hits, selectBuf);
//}
//void Mouse(int button, int state, int cursorX, int cursorY)
//{
//    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
//    startPicking(cursorX, cursorY);
//    //just check if the cursor is in the circle...
//}
//int main (int argc,  char *argv[])
//{
//    glutInit(&argc, argv);
//    glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
//    glutInitWindowSize(480,480);
//    int windowHandle = glutCreateWindow("Force layout");
//    
//    glutDisplayFunc(redraw);
//    glutReshapeFunc(reshape);
//    glutKeyboardFunc(key);
//    glutMouseFunc(Mouse);
//    glutIdleFunc(idle);
//    
//    CGLContextObj ctx = CGLGetCurrentContext();
//    const int interval = 0;
//    
//    CGLSetParameter(ctx, kCGLCPSwapInterval, &interval);
//    
//    glutMainLoop();
//    return 0;
//}


#include <stdio.h>
#include <GLUT/GLUT.h>
#include <math.h>

#define DISTANCE(X,Y,CX,CY) sqrt((X-CX)*(X-CX)+(Y-CY)*(Y-CY))

const int n = 1000;
const GLfloat r = 0.5f;
const GLfloat cx = 1.0f, cy = 1.0f;

const GLint windowX = 600, windowY = 600;
const GLfloat scaleTime = 5;
bool selected = false;

void coorTrans(const int wx, const int wy, float& x, float& y)
{
    x = (wx-windowX/2)/(windowX/2/scaleTime);
    y = -(wy-windowY/2)/(windowY/2/scaleTime);
}
void DrawCircle()
{
    int i;
//    glBegin(GL_LINE_LOOP);
    glPushMatrix();
    if(!selected)
        glColor3f(1.0, 1.0, 1.0);
    else
        glColor3f(1.0, 1.0, 0.0);
    glPopMatrix();
    glBegin(GL_POLYGON);
    for(i=0; i<n; ++i)
        glVertex2f(cx+r*cos(2*M_PI/n*i), cy+r*sin(2*M_PI/n*i));
    glEnd();
    glFlush();
}
void Mouse(int button, int state, int cursorX, int cursorY)
{
    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
    {
        printf("%d, %d\n", cursorX, cursorY);
        float x = 0, y = 0;
        coorTrans(cursorX, cursorY, x, y);
        if(DISTANCE(x, y, cx, cy) < r)
        {
            selected = !selected;
            DrawCircle();
        }
    }
}
void myDisplay()
{
    glClear(GL_COLOR_BUFFER_BIT);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glScalef(1.0/scaleTime, 1.0/scaleTime, 1.0/scaleTime);
    DrawCircle();
}
void initFunc()
{
    
}

int main(int argc, char *argv[])
{
    glutInit(&argc, argv);
    initFunc();
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowPosition(100, 100);
    glutInitWindowSize(windowX, windowY);
    glutCreateWindow("Draw a circle");
    glutDisplayFunc(myDisplay);
    glutMouseFunc(Mouse);
    glutMainLoop();
    return 0;
}
