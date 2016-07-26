#ifndef WINDOW_H
#define WINDOW_H
#include <QOpenGLWidget>
#include <glu.h>
#include <OpenGL.h>

#include <boost/graph/fruchterman_reingold.hpp>
#include <boost/graph/random_layout.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topology.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/progress.hpp>
#include <boost/shared_ptr.hpp>
#include <stdio.h>
#include <GLUT/GLUT.h>
#include <OpenGL/OpenGL.h>
#include <GLUT/gutil.h>
#include <OpenGL/glext.h>
#include <utility>
#include <vector>
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <algorithm>
#include <QKeyEvent>
#include <QMouseEvent>

#include "betweenness.h"
#include "graph.h"
#include "define.h"

typedef boost::rectangle_topology<> topology_type;
typedef topology_type::point_type point_type;

typedef adjacency_list<listS, vecS, undirectedS,
property<vertex_name_t, std::string> > Graph;

typedef graph_traits<Graph>::vertex_descriptor Vertex;

typedef std::map<std::string, Vertex> NameToVertex;

class Window : public QOpenGLWidget {
  Q_OBJECT
public:
    Window();

    virtual void initializeGL();
    virtual void paintGL();
    virtual void resizeGL(int width, int height);
    virtual void keyPressEvent(QKeyEvent *event);
    virtual void mousePressEvent(QMouseEvent *event);
//    virtual void mouseMoveEvent(QMouseEvent *event);

    void updateWindow();
    void updateData();
    void updateColor();

protected:
    bool datachanged, colorchanged;

//variables added by ansjcy
private:
    const int circlePoints = 100;
    const GLfloat r = 0.008f;
    GLint nodeSize = 6;
    std::vector<bool> selected;
    int selectNode = -1;
    GLint edgeSize = 5;
    std::vector<bool> selectedEdge;
    std::vector<std::pair<GLfloat, GLfloat> > myNodes;
    std::vector<std::pair<GLint, GLint> > myEdges;
    Betweenness bc;
    std::vector<Node*> nodes;
    std::vector<Edge*> edges;
    const GLint windowX = 850, windowY = 850;
    GLfloat scaleTime = 1;
    GLfloat scaleTotal = 1;
    GLfloat transformX = 0;
    GLfloat transformY = 0;
    GLfloat transformXTotal = 0;
    GLfloat transformYTotal = 0;
    bool scaleFlag = false;
    bool transformFlag = false;



//functions added by ansjcy
private:
    Vertex get_vertex(const std::string& name, Graph& g, NameToVertex& names);
    void initFunc();
    void coorTrans(const int wx, const int wy, float& x, float& y);
    void DrawCircle(float cx, float cy, float r, int num_segments);
    std::vector<std::string> getNextLineAndSplitIntoTokens(std::istream& str);

signals:
    void emitUpdate();
    void emitColorUpdate();
    void emitDataUpdate();

};

#endif // WINDOW_H
