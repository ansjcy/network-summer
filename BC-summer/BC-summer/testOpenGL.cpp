//#include <stdio.h>
//#include <iostream>
//#include <GLUT/GLUT.h>
//#include <utility>
//#include <vector>
//#include <math.h>
//
//#define DISTANCE(X,Y,CX,CY) sqrt((X-CX)*(X-CX)+(Y-CY)*(Y-CY))
//
//const int circlePoints = 1000;
//const GLfloat r = 0.5f;
//
//GLint nodeSize = 6;
//std::vector<bool> selected;
//GLint edgeSize = 5;
//std::vector<bool> selectedEdge;
//
//
//std::vector<std::pair<GLfloat, GLfloat> > myNodes;
//std::vector<std::pair<GLint, GLint> > myEdges;
//
//
//const GLint windowX = 600, windowY = 600;
//const GLfloat scaleTime = 5;
//
//
//void coorTrans(const int wx, const int wy, float& x, float& y)
//{
//    x = (wx-windowX/2)/(windowX/2/scaleTime);
//    y = -(wy-windowY/2)/(windowY/2/scaleTime);
//}
//
//
//void DrawCircle()
//{
//    for(int i = 0; i < myNodes.size(); i++)
//    {
//        glPushMatrix();
//        if(!selected[i])
//            glColor3f(1.0, 1.0, 1.0);
//        else
//            glColor3f(1.0, 1.0, 0.0);
//        glPopMatrix();
//        glBegin(GL_POLYGON);
//        for(int k=0; k<circlePoints; k++)
//            glVertex2f(myNodes[i].first+r*cos(2*M_PI/circlePoints*k), myNodes[i].second+r*sin(2*M_PI/circlePoints*k));
//        glEnd();
//    }
//    
//    for(int i = 0; i < myEdges.size(); i++)
//    {
//        glPushMatrix();
//        if(!selectedEdge[i])
//            glColor3f(1.0, 1.0, 1.0);
//        else
//            glColor3f(1.0, 1.0, 0.0);
//        glPopMatrix();
//        
//        glBegin(GL_LINE_STRIP);
//        GLfloat startX = myNodes[myEdges[i].first].first, startY = myNodes[myEdges[i].first].second;
//        GLfloat endX = myNodes[myEdges[i].second].first, endY = myNodes[myEdges[i].second].second;
//        
//        glVertex2f(startX,startY);
//        glVertex2f(endX,endY);
//        glEnd();
//    }
//    glFlush();
//}
//
//
//void Mouse(int button, int state, int cursorX, int cursorY)
//{
//    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
//    {
//        printf("%d, %d\n", cursorX, cursorY);
//        float x = 0, y = 0;
//        for(int i = 0; i < myNodes.size(); i++)
//        {
//            coorTrans(cursorX, cursorY, x, y);
//            if(DISTANCE(x, y, myNodes[i].first, myNodes[i].second) < r)
//            {
//                printf("i = %d\n", i);
//                selected[i] = !selected[i];
//                DrawCircle();
//            }
//        }
//        
//    }
//}
//void myDisplay()
//{
//    glClear(GL_COLOR_BUFFER_BIT);
//    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
//    glScalef(1.0/scaleTime, 1.0/scaleTime, 1.0/scaleTime);
//    DrawCircle();
//}
//void initFunc()
//{
//    myNodes.push_back(std::make_pair(0.0, 0.0));
//    myNodes.push_back(std::make_pair(-3.0, -2.5));
//    myNodes.push_back(std::make_pair(0.0, -4.0));
//    myNodes.push_back(std::make_pair(2.0, 2.0));
//    myNodes.push_back(std::make_pair(1.5, 4.0));
//    myNodes.push_back(std::make_pair(3.5, 0.0));
//    
//    myEdges.push_back(std::make_pair(0, 1));
//    myEdges.push_back(std::make_pair(0, 2));
//    myEdges.push_back(std::make_pair(0, 3));
//    myEdges.push_back(std::make_pair(0, 4));
//    myEdges.push_back(std::make_pair(3, 5));
//    
//    for(int i = 0; i < myNodes.size(); i++)
//    {
//        selected.push_back(false);
//    }
//    for(int i = 0; i < myEdges.size(); i++)
//    {
//        selectedEdge.push_back(false);
//    }
//}
//
//int main(int argc, char *argv[])
//{
//    glutInit(&argc, argv);
//    initFunc();
//    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH);
//    glutInitWindowPosition(100, 100);
//    glutInitWindowSize(windowX, windowY);
//    glutCreateWindow("Draw a circle");
//    glutDisplayFunc(myDisplay);
//    glutMouseFunc(Mouse);
//    glutMainLoop();
//    return 0;
//}