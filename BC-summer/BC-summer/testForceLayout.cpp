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
#include <utility>
#include <vector>
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <algorithm>

#include "betweenness.hpp"
#include "graph.hpp"
//

using namespace boost;

typedef boost::rectangle_topology<> topology_type;
typedef topology_type::point_type point_type;

typedef adjacency_list<listS, vecS, undirectedS,
property<vertex_name_t, std::string> > Graph;

typedef graph_traits<Graph>::vertex_descriptor Vertex;

typedef std::map<std::string, Vertex> NameToVertex;

Vertex get_vertex(const std::string& name, Graph& g, NameToVertex& names)
{
    NameToVertex::iterator i = names.find(name);
    if (i == names.end())
        i = names.insert(std::make_pair(name, add_vertex(name, g))).first;
    return i->second;
}

class progress_cooling : public linear_cooling<double>
{
    typedef linear_cooling<double> inherited;
    
public:
    explicit progress_cooling(std::size_t iterations) : inherited(iterations)
    {
        display.reset(new progress_display(iterations + 1, std::cerr));
    }
    
    double operator()()
    {
        ++(*display);
        return inherited::operator()();
    }
    
private:
    shared_ptr<boost::progress_display> display;
};
//******* for openGL functions *************



#define DISTANCE(X,Y,CX,CY) sqrt((X-CX)*(X-CX)+(Y-CY)*(Y-CY))

const int circlePoints = 100;
const GLfloat r = 0.01f;

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


const GLint windowX = 800, windowY = 800;
const GLfloat scaleTime = 1;



void coorTrans(const int wx, const int wy, float& x, float& y)
{
    x = (wx-windowX/2)/(windowX/2/scaleTime);
    y = -(wy-windowY/2)/(windowY/2/scaleTime);
}

//double colorStart[] = {1.0, 0.0, 0.0};
//double colorEnd[] = {0.0, 0.0, 1.0};

void DrawCircle()
{
    
    for(int i = 0; i < myEdges.size(); i++)
    {
        glPushMatrix();
        
        if(selectNode < 0)
        {
//            glColor4f(0.0, 0.0, 0.0, 1.0);
            auto minMaxvalue = std::minmax_element(std::begin(nodes[myEdges[i].first]->sensitivityValues), std::end(nodes[myEdges[i].first]->sensitivityValues));
            double red = 0.0, blue = 0.0;
            double sensitivity = nodes[myEdges[i].first]->sensitivityValues[myEdges[i].second];
            if(sensitivity < 0)
            {
                //red
                red = sensitivity/(*minMaxvalue.first)+0.2;
            }
            else
            {
                blue = sensitivity/(*minMaxvalue.second)+0.2;
            }
            glColor4f(red, 0.0, blue, 1.0);
            
        }
        else
            glColor4f(0.1, 0.1, 0.1, 1.0);
        
        
        
        //color the edges according to the sensitivity value
        if(selectNode == myEdges[i].first)
        {
            auto minMaxvalue = std::minmax_element(std::begin(nodes[selectNode]->sensitivityValues), std::end(nodes[selectNode]->sensitivityValues));
            double red = 0.0, blue = 0.0;
            double sensitivity = nodes[selectNode]->sensitivityValues[myEdges[i].second];
            
            if(sensitivity < 0)
            {
                //red
                red = sensitivity/(*minMaxvalue.first)+0.2;
            }
            else
            {
                blue = sensitivity/(*minMaxvalue.second)+0.2;
            }
            
            glColor4f(red, 0.0, blue, 1.0);
        }
        
        glPopMatrix();
        
        glBegin(GL_LINE_STRIP);
        GLfloat startX = myNodes[myEdges[i].first].first, startY = myNodes[myEdges[i].first].second;
        GLfloat endX = myNodes[myEdges[i].second].first, endY = myNodes[myEdges[i].second].second;
        
        glVertex2f(startX,startY);
        glVertex2f(endX,endY);
        glEnd();
    }
    
    
    for(int i = 0; i < myNodes.size(); i++)
    {
        glPushMatrix();
//        if(!selected[i])
//            glColor4f(1.0, 1.0, 1.0, 1.0);
//        else
//            glColor4f(1.0, 1.0, 0.0, 1.0);
        if(selectNode <= 0)
            glColor4f(0.1, 0.1, 0.1, 1.0);
        
        else
        {
            double signValue = nodes[selectNode]->sensitivityValues[i];
            if(signValue < 0)
                glColor4f(1.0, 0.0, 0.0, 1.0);
            else
                glColor4f(0.0, 0.0, 1.0, 1.0);
            
        }
        glPopMatrix();
        //
        glBegin(GL_POLYGON);
        
        
        
        //gldrawarray.. gltranangle
        for(int k=0; k<circlePoints; k++)
            glVertex2f(myNodes[i].first+r*cos(2*M_PI/circlePoints*k), myNodes[i].second+r*sin(2*M_PI/circlePoints*k));
        glEnd();
    }
    
    
    glFlush();
}


void Mouse(int button, int state, int cursorX, int cursorY)
{
    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
    {
        printf("%d, %d\n", cursorX, cursorY);
        float x = 0, y = 0;
        for(int i = 0; i < myNodes.size(); i++)
        {
            coorTrans(cursorX, cursorY, x, y);
            if(DISTANCE(x, y, myNodes[i].first, myNodes[i].second) < r)
            {
                printf("i = %d\n", i);
                selected[i] = !selected[i];
                if(selectNode == i)
                    selectNode = -1;
                else
                    selectNode = i;
                DrawCircle();
            }
        }
        
    }
}
void myDisplay()
{
    glClear(GL_COLOR_BUFFER_BIT);
    glClearColor(0.9f, 0.9f, 0.9f, 0.0f);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glEnable(GL_BLEND);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glScalef(1.0/scaleTime, 1.0/scaleTime, 1.0/scaleTime);
    DrawCircle();
}

std::vector<std::string> getNextLineAndSplitIntoTokens(std::istream& str)
{
    std::vector<std::string>   result;
    std::string                line;
    std::getline(str,line);
    
    std::stringstream          lineStream(line);
    std::string                cell;
    
    while(std::getline(lineStream,cell, ','))
    {
        result.push_back(cell);
    }
    return result;
}

void initFunc()
{

    std::ifstream nodeFile, edgeFile;
    nodeFile.open("/Users/anakin/Downloads/data/adjnoun.nodes.csv");
    edgeFile.open("/Users/anakin/Downloads/data/adjnoun.edges.csv");
//    nodeFile.open("/Users/anakin/Downloads/data/celegansneural.nodes.csv");
//    edgeFile.open("/Users/anakin/Downloads/data/celegansneural.edges.csv");
    
    std::vector<std::string> result = getNextLineAndSplitIntoTokens(nodeFile);
    while ((result = getNextLineAndSplitIntoTokens(nodeFile)).size()!=0) {
        myNodes.push_back(std::make_pair(std::stod(result[1]), std::stod(result[2])));
    }
    
    
    result = getNextLineAndSplitIntoTokens(edgeFile);
    while ((result = getNextLineAndSplitIntoTokens(edgeFile)).size()!=0) {
        myEdges.push_back(std::make_pair(std::stod(result[0]), std::stod(result[1])));
    }
    
    for(int i = 0; i < myNodes.size(); i++)
    {
        selected.push_back(false);
    }
    for(int i = 0; i < myEdges.size(); i++)
    {
        selectedEdge.push_back(false);
    }
}




void idle()
{
    glutPostRedisplay();
}
int main(int argc, char* argv[])
{

    initFunc();
    
    
//*************** do the force layout **************************
    
    int iterations = 100;
    
    double width = 2.0;
    double height = 2.0;
    Graph g;
    NameToVertex names;
    
    for(int i = 0; i < myEdges.size(); i++)
    {
        add_edge(get_vertex(std::to_string(myEdges[i].first), g, names), get_vertex(std::to_string(myEdges[i].second), g, names), g);
        
    }
    
    typedef std::vector<point_type> PositionVec;
    PositionVec position_vec(num_vertices(g));
    typedef iterator_property_map<PositionVec::iterator,
    property_map<Graph, vertex_index_t>::type>
    PositionMap;
    PositionMap position(position_vec.begin(), get(vertex_index, g));
    
    minstd_rand gen;
    topology_type topo(gen, -width/2, -height/2, width/2, height/2);
    random_graph_layout(g, position, topo);
    fruchterman_reingold_force_directed_layout
    (g, position, topo,
     cooling(progress_cooling(iterations)));
    
    graph_traits<Graph>::vertex_iterator vi, vi_end;
    for (boost::tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi) {
        std::cout << get(vertex_name, g, *vi) << '\t'
        << position[*vi][0] << '\t' << position[*vi][1] << std::endl;
    }
    for(boost::tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
    {
        myNodes[std::stoi(get(vertex_name, g, *vi))].first = position[*vi][0];
        myNodes[std::stoi(get(vertex_name, g, *vi))].second = position[*vi][1];
    }
    
    
//******** for bc ***********

    for(int i = 0; i < myNodes.size(); i++)
    {
        Node* tmpNode = new Node();
        tmpNode->setIndex(i);
        nodes.push_back(tmpNode);
    }
    for(int i = 0; i < myEdges.size(); i++)
    {
        nodes[myEdges[i].first]->addEdge(nodes[myEdges[i].second]);
        nodes[myEdges[i].second]->addEdge(nodes[myEdges[i].first]);
    }

    bc.compute(nodes, true);

//***************************


    //for openGL functions
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowPosition(100, 100);
    glutInitWindowSize(windowX, windowY);
    glutCreateWindow("Draw a circle");
    glutDisplayFunc(myDisplay);
    glutMouseFunc(Mouse);
    glutIdleFunc(idle);
    glutMainLoop();
    

    
    
    
    
    return 0;
}