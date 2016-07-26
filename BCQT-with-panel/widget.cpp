#include "widget.h"
#include "ui_widget.h"

using namespace boost;

enum Mode{
    CENTRALITY,
    SENSITIVITY_MEAN,
    SENSITIVITY_VARIANCE
} MODE = CENTRALITY;


Widget::Widget(QOpenGLWidget *parent) :
    QOpenGLWidget(parent),
    ui(new Ui::Widget)
{
    ui->setupUi(this);

    QSurfaceFormat format;
    format.setSamples(8);
    format.setStencilBufferSize(8);
    setFormat(format);
    printf("Version Major:%d minor:%d \n",format.version().first, format.version().second);
    datachanged = colorchanged = true;
}

Widget::~Widget()
{
    delete ui;
}

void Widget::updateWindow(){
  datachanged = colorchanged = true;
  update();
}

void Widget::updateData(){
  datachanged = true;
  update();
}

void Widget::updateColor(){
  colorchanged = true;
  update();
}

//*********** helper functions *****************
std::vector<std::string> Widget::getNextLineAndSplitIntoTokens(std::istream& str)
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

Vertex Widget::get_vertex(const std::string& name, Graph& g, NameToVertex& names)
{
    NameToVertex::iterator i = names.find(name);
    if (i == names.end())
        i = names.insert(std::make_pair(name, add_vertex(name, g))).first;
    return i->second;
}

/*
 * transfom window such as 0 - 800 to -1 -> 1 (divided by the scale factor)
 *      |y
 *      |
 * ---------->x
 *      |
 *      |
 */
void Widget::coorTrans(const int wx, const int wy, float& x, float& y)
{
    x = (wx-windowX/2)/(windowX/2/scaleTotal) - transformXTotal;
    y = -((wy-windowY/2)/(windowY/2/scaleTotal) + transformYTotal);
}

//void Widget::DrawCircle(float cx, float cy, float r, int num_segments)
//{
//    float theta = 2 * M_PI / float(num_segments);
//    float c = cosf(theta);//precalculate the sine and cosine
//    float s = sinf(theta);
//    float t;

//    float x = r;//we start at angle = 0
//    float y = 0;

//    glBegin(GL_LINE_LOOP);
//    for(int ii = 0; ii < num_segments; ii++)
//    {
//        glVertex2f(x + cx, y + cy);//output vertex
//        //apply the rotation matrix
//        t = x;
//        x = c * x - s * y;
//        y = s * t + c * y;
//    }
//    glEnd();
//}

void Widget::initFunc()
{

    this->setFocus();

    std::ifstream nodeFile, edgeFile;
    nodeFile.open("/Users/anakin/Downloads/data/adjnoun.nodes.csv");
    edgeFile.open("/Users/anakin/Downloads/data/adjnoun.edges.csv");

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
void Widget::keyPressEvent(QKeyEvent *event)
{
    scaleFlag = false;
    transformFlag = false;
    std::cout << "key::" << event->key() << std::endl;
    switch (event->key())
        {
    case Qt::Key_Q:
        MODE = CENTRALITY;
        break;
    case Qt::Key_W:
        std::cout << "hereW" << std::endl;
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
    case Qt::Key_Escape:
        this->close();
            default:
                break;
        }
    std::cout << "before repaint:: " << scaleTime << std::endl;
//    this->repaint();
    this->update();

}
void Widget::mousePressEvent(QMouseEvent *event)
{
    this->setFocus();
    if(event->button() == Qt::LeftButton)
    {
         std::cout << event->x() << " " << event->y() << std::endl;
         float x = 0, y = 0;
         for(int i = 0; i < myNodes.size(); i++)
         {
             coorTrans(event->x(), event->y(), x, y);
             if(DISTANCE(x, y, myNodes[i].first, myNodes[i].second) < r*nodeSize)
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
//         this->repaint();
         this->update();
     }

    QOpenGLWidget::mousePressEvent(event);
}

void Widget::resizeGL(int width, int height){
  (void) width;
  (void) height;
}

void Widget::initializeGL(){
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
void Widget::drawACircle(float cx, float cy, float r, rgb colors, float alpha)
{
    int triangleAmount = 20; //# of triangles used to draw circle

        //GLfloat radius = 0.8f; //radius
    GLfloat twicePi = 2.0f * M_PI;
    GLfloat x = cx, y = cy;

    glColor4f(colors.r/255, colors.g/255, colors.b/255, alpha);

    glBegin(GL_TRIANGLE_FAN);
        glVertex2f(x, y); // center of circle
        for(int k = 0; k <= triangleAmount;k++) {
            glVertex2f(
                    x + (r * nodeSize * cos(k *  twicePi / triangleAmount)),
                y + (r * nodeSize * sin(k * twicePi / triangleAmount))
            );
        }
    glEnd();
}
void Widget::drawALine(float startX, float startY, float endX, float endY, float lineWidth, rgb colors, float alpha)
{
    glLineWidth(lineWidth);
    glColor4f(colors.r/255, colors.g/255, colors.b/255, alpha);
    glBegin(GL_LINE_STRIP);

    glVertex2f(startX,startY);
    glVertex2f(endX,endY);
    glEnd();
}

void Widget::paintGL(){
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

    rgb whiteValue;
    whiteValue.r = 229;
    whiteValue.g = 229;
    whiteValue.b = 229;

    //draw edge with sensitivity = 0 first
    for(int i = 0; i < myEdges.size(); i++)
    {
        if(fabs(nodes[myEdges[i].first]->sensitivityValues[myEdges[i].second]-0) < 1e-9 || (selectNode >= 0 && selectNode != myEdges[i].first))
        {
            drawALine(myNodes[myEdges[i].first].first, myNodes[myEdges[i].first].second, myNodes[myEdges[i].second].first, myNodes[myEdges[i].second].second,
                    edgeSize, whiteValue, edgeTran);
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

            drawALine(myNodes[myEdges[i].first].first, myNodes[myEdges[i].first].second, myNodes[myEdges[i].second].first, myNodes[myEdges[i].second].second,
                    edgeSize, colors, edgeTran);
        }
    }

    //draw shadows
    for(int i = 0; i < myNodes.size(); i++)
    {
        rgb colors;
        colors.r = 0.3;
        colors.g = 0.3;
        colors.b = 0.3;

        drawACircle(myNodes[i].first+r*0.2, myNodes[i].second-r*0.2, r*nodeSize, colors, 0.3*nodeTran);
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
        rgb colors;
        if(selectNode < 0)
        {

            double value = centralityValues[i];
            if(MODE == SENSITIVITY_MEAN)
                value = sensitivityMeans[i];
            else if(MODE == SENSITIVITY_VARIANCE)
                value = sensitivityVariances[i];

            colors = getRGBValue(*minMaxvalue.first, *minMaxvalue.second, value);
            //glColor4f(colors.r/255, colors.g/255, colors.b/255, 1);

        }
        else
        {
            auto minMaxvalue = std::minmax_element(std::begin(nodes[selectNode]->sensitivityValues), std::end(nodes[selectNode]->sensitivityValues));
            double sensitivity = nodes[selectNode]->sensitivityValues[i];

            colors = getRGBValue(*minMaxvalue.first, *minMaxvalue.second, sensitivity);
            //glColor4f(colors.r/255, colors.g/255, colors.b/255, 1);

        }

        drawACircle(myNodes[i].first, myNodes[i].second, r*nodeSize, colors, nodeTran);

    }
    glFlush();
}


void Widget::on_nodeTran_valueChanged(double arg1)
{
    nodeTran = arg1;
    this->update();
}

void Widget::on_nodeSize_valueChanged(double arg1)
{
    nodeSize = arg1;
    this->update();
}

void Widget::on_edgeTran_valueChanged(double arg1)
{
    edgeTran = arg1;
    this->update();
}

void Widget::on_edgeSize_valueChanged(double arg1)
{
    edgeSize = arg1;
    this->update();
}
