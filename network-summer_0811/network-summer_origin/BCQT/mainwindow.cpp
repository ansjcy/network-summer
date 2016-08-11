#include "mainwindow.h"


using namespace boost;

enum Mode{
    CENTRALITY,
    SENSITIVITY_MEAN,
    SENSITIVITY_VARIANCE
} MODE = CENTRALITY;



Window::Window()
  : QOpenGLWidget()
{
  QSurfaceFormat format;
  format.setSamples(8);
  format.setStencilBufferSize(8);
  setFormat(format);
  printf("Version Major:%d minor:%d \n",format.version().first, format.version().second);
  datachanged = colorchanged = true;
}

void Window::updateWindow(){
  datachanged = colorchanged = true;
  update();
}

void Window::updateData(){
  datachanged = true;
  update();
}

void Window::updateColor(){
  colorchanged = true;
  update();
}

//*********** helper functions *****************
std::vector<std::string> Window::getNextLineAndSplitIntoTokens(std::istream& str)
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

Vertex Window::get_vertex(const std::string& name, Graph& g, NameToVertex& names)
{
    NameToVertex::iterator i = names.find(name);
    if (i == names.end())
        i = names.insert(std::make_pair(name, add_vertex(name, g))).first;
    return i->second;
}

void Window::coorTrans(const int wx, const int wy, float& x, float& y)
{
    x = (wx-windowX/2)/(windowX/2/scaleTotal);
    y = -(wy-windowY/2)/(windowY/2/scaleTotal);
}

void Window::DrawCircle(float cx, float cy, float r, int num_segments)
{
    float theta = 2 * M_PI / float(num_segments);
    float c = cosf(theta);//precalculate the sine and cosine
    float s = sinf(theta);
    float t;

    float x = r;//we start at angle = 0
    float y = 0;

    glBegin(GL_LINE_LOOP);
    for(int ii = 0; ii < num_segments; ii++)
    {
        glVertex2f(x + cx, y + cy);//output vertex
        //apply the rotation matrix
        t = x;
        x = c * x - s * y;
        y = s * t + c * y;
    }
    glEnd();
}

void Window::initFunc()
{

    std::ifstream nodeFile, edgeFile;
    nodeFile.open("/Users/anakin/Downloads/data/serengeti-foodweb.nodes.csv");
    edgeFile.open("/Users/anakin/Downloads/data/serengeti-foodweb.edges.csv");

//    nodeFile.open("/Users/anakin/Downloads/data/netscience.nodes.csv");
//    edgeFile.open("/Users/anakin/Downloads/data/netscience.edges.csv");

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


//*************** main functions *******************
void Window::keyPressEvent(QKeyEvent *event)
{
    scaleFlag = false;
    transformFlag = false;
    switch (event->key())
        {
    case Qt::Key_Q:
        MODE = CENTRALITY;
        break;
    case Qt::Key_W:
        MODE = SENSITIVITY_MEAN;
        break;
    case Qt::Key_E:
        MODE = SENSITIVITY_VARIANCE;
        break;
    case Qt::Key_U:
        scaleFlag = true;
        scaleTime = 1.1;
        break;
    case Qt::Key_O:
        scaleFlag =true;
        scaleTime = 0.9;
        break;
    case Qt::Key_I:
        transformFlag = true;
        transformX = 0;
        transformY = 0.1;
        break;
    case Qt::Key_K:
        transformFlag = true;
        transformX = 0;
        transformY = -0.1;
        break;
    case Qt::Key_J:
        transformFlag = true;
        transformX = -0.1;
        transformY = 0;
        break;
    case Qt::Key_L:
        transformFlag = true;
        transformX = 0.1;
        transformY = 0;
        break;
            default:
                break;
        }
    std::cout << "before repaint:: " << scaleTime << std::endl;
//    this->repaint();
    this->update();

}
void Window::mousePressEvent(QMouseEvent *event)
{
    if(event->button() == Qt::LeftButton)
    {
         std::cout << event->x() << " " << event->y() << std::endl;
         float x = 0, y = 0;
         for(int i = 0; i < myNodes.size(); i++)
         {
             coorTrans(event->x(), event->y(), x, y);
             if(DISTANCE(x, y, myNodes[i].first, myNodes[i].second) < r)
             {
                 printf("i = %d\n", i);
                 selected[i] = !selected[i];
                 if(selectNode == i)
                     selectNode = -1;
                 else
                     selectNode = i;
                 std::cout << "mean: " << nodes[i]->sensitivityMean << std::endl;
             }
         }
         transformFlag = false;
         scaleFlag = false;
         this->repaint();
//         this->update();
     }

    QOpenGLWidget::mousePressEvent(event);
}

void Window::resizeGL(int width, int height){
  (void) width;
  (void) height;
}

void Window::initializeGL(){
    glClear(GL_COLOR_BUFFER_BIT);
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glEnable(GL_BLEND);
    //smooth
    glEnable(GL_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


    glClearDepthf(1.0f);
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_DEPTH_TEST);

    initFunc();

    //******* do the force layout *******
    int iterations = 1000;

    double width = 1.8;
    double height = 1.8;
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
    (g, position, topo);

    graph_traits<Graph>::vertex_iterator vi, vi_end;
#ifdef DEBUG_FORCE_LAYOUT_POS
    for (boost::tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi) {
        std::cout << get(vertex_name, g, *vi) << '\t'
        << position[*vi][0] << '\t' << position[*vi][1] << std::endl;
    }
#endif
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
}

void Window::paintGL(){
    std::cout << "scaleTime::" << scaleTime << std::endl;

    if(scaleFlag)
    {
        glScalef(1.0/scaleTime, 1.0/scaleTime, 1.0/scaleTime);
        scaleTotal *= scaleTime;
    }
    if(transformFlag)
    {
        glTranslatef(transformX, transformY,0.0f);
        transformXTotal += transformX;
        transformYTotal += transformY;
    }

    double whiteValue = 0.9;

    //draw edge with sensitivity = 0 first
    for(int i = 0; i < myEdges.size(); i++)
    {
        if(fabs(nodes[myEdges[i].first]->sensitivityValues[myEdges[i].second]-0) < 1e-9 || (selectNode >= 0 && selectNode != myEdges[i].first))
        {
            glColor4f(whiteValue, whiteValue, whiteValue, 1);
            glBegin(GL_LINE_STRIP);
            GLfloat startX = myNodes[myEdges[i].first].first, startY = myNodes[myEdges[i].first].second;
            GLfloat endX = myNodes[myEdges[i].second].first, endY = myNodes[myEdges[i].second].second;
            glVertex2f(startX,startY);
            glVertex2f(endX,endY);
            glEnd();

        }
    }

    for(int i = 0; i < myEdges.size(); i++)
    {
        if(fabs(nodes[myEdges[i].first]->sensitivityValues[myEdges[i].second]-0) < 1e-9)
            continue;

        if(selectNode < 0 || selectNode == myEdges[i].first)
        {
            auto minMaxvalue = std::minmax_element(std::begin(nodes[myEdges[i].first]->sensitivityValues), std::end(nodes[myEdges[i].first]->sensitivityValues));
            double sensitivity = nodes[myEdges[i].first]->sensitivityValues[myEdges[i].second];

            rgb colors = getRGBValue(*minMaxvalue.first, *minMaxvalue.second, sensitivity);
            glColor4f(colors.r/255, colors.g/255, colors.b/255, 0.6);


            glBegin(GL_LINE_STRIP);
            GLfloat startX = myNodes[myEdges[i].first].first, startY = myNodes[myEdges[i].first].second;
            GLfloat endX = myNodes[myEdges[i].second].first, endY = myNodes[myEdges[i].second].second;

            glVertex2f(startX,startY);
            glVertex2f(endX,endY);
            glEnd();

        }
    }

    //draw shadows
    for(int i = 0; i < myNodes.size(); i++)
    {
        glColor4f(0.3, 0.3, 0.3, 0.3);
        glBegin(GL_POLYGON);

        //gldrawarray.. gltranangle
        for(int k=0; k<circlePoints; k++)
            glVertex2f(myNodes[i].first+r*0.2+r*cos(2*M_PI/circlePoints*k), myNodes[i].second-r*0.2+r*sin(2*M_PI/circlePoints*k));
        glEnd();

        DrawCircle(myNodes[i].first+r*0.2, myNodes[i].second-r*0.2, r*1.1, circlePoints);
    }

    std::vector<double> sensitivityMeans;
    std::vector<double> sensitivityVariances;
    std::vector<double> centralityValues;
    for(int i = 0; i < nodes.size(); i++)
    {
        sensitivityMeans.push_back(nodes[i]->sensitivityMean);
        sensitivityVariances.push_back(nodes[i]->sensitivityVariance);
        centralityValues.push_back(nodes[i]->centralityValue);
    }

    //default mode is centrality mode
    auto minMaxvalue = std::minmax_element(std::begin(centralityValues), std::end(centralityValues));
    if(MODE == SENSITIVITY_MEAN)
        minMaxvalue = std::minmax_element(std::begin(sensitivityMeans), std::end(sensitivityMeans));
    else if(MODE == SENSITIVITY_VARIANCE)
        minMaxvalue = std::minmax_element(std::begin(sensitivityVariances), std::end(sensitivityVariances));

    for(int i = 0; i < myNodes.size(); i++)
    {
        if(selectNode < 0)
        {

            double value = centralityValues[i];
            if(MODE == SENSITIVITY_MEAN)
                value = sensitivityMeans[i];
            else if(MODE == SENSITIVITY_VARIANCE)
                value = sensitivityVariances[i];

            rgb test = getRGBValue(*minMaxvalue.first, *minMaxvalue.second, value);
            glColor4f(test.r/255, test.g/255, test.b/255, 1);

        }
        else
        {
            auto minMaxvalue = std::minmax_element(std::begin(nodes[selectNode]->sensitivityValues), std::end(nodes[selectNode]->sensitivityValues));
            double sensitivity = nodes[selectNode]->sensitivityValues[i];

            rgb colors = getRGBValue(*minMaxvalue.first, *minMaxvalue.second, sensitivity);
            glColor4f(colors.r/255, colors.g/255, colors.b/255, 1);

        }
        //
//        glBegin(GL_POLYGON);

//        //gldrawarray.. gltranangle
//        for(int k=0; k<circlePoints; k++)
//            glVertex2f(myNodes[i].first+r*cos(2*M_PI/circlePoints*k), myNodes[i].second+r*sin(2*M_PI/circlePoints*k));
//        glEnd();

//        DrawCircle(myNodes[i].first, myNodes[i].second, r*1.1, circlePoints);

        int k;
        int triangleAmount = 20; //# of triangles used to draw circle

            //GLfloat radius = 0.8f; //radius
        GLfloat twicePi = 2.0f * M_PI;
        GLfloat x = myNodes[i].first, y = myNodes[i].second;

            glBegin(GL_TRIANGLE_FAN);
                glVertex2f(x, y); // center of circle
                for(k = 0; k <= triangleAmount;k++) {
                    glVertex2f(
                            x + (r * cos(k *  twicePi / triangleAmount)),
                        y + (r * sin(k * twicePi / triangleAmount))
                    );
                }
            glEnd();

    }
    glFlush();
}
