//
//  graph.hpp
//  BC
//
//  Created by Anakin on 16/7/14.
//  Copyright © 2016年 Anakin. All rights reserved.
//

#ifndef graph_hpp
#define graph_hpp

#include <string>
#include <vector>
#include <map>
#include <cmath>
//#include "properties.hpp"
class Edge;
class Node
{
public:
    std::vector<Edge*> edges;
    std::map<Node*, double> pred;
    int index;
    double centralityValue;
    std::vector<double> sensitivityValues;
    double sensitivityMean;
    double sensitivityVariance;
//    double *sensitivityValues;
    
//    std::string property;
    Node(){ index = -1; }
    ~Node()
    {
        for(int i = 0; i < edges.size(); i++)
            delete edges[i];
    }
    
    Edge* addEdge(Node* n, double weight);
    int findEdge(Node* n);
    Edge* getEdge(int idx);
    void setIndex(int n);
    int getIndex();
    bool operator==(Node &b);
    bool operator==(int idx);
    
};
class Edge
{
public:
    int multiplicity;
    float weight;
    double centralityValue;
    Node* node0;
    Node* node1;
    
    Edge(){ node0 = 0; node1 = 0; multiplicity = 0; weight = 0; }
    void setNode0(Node* n){ node0 = n; }
    void setNode1(Node* n){ node1 = n; }
    Node* getNode0(){return node0;}
    Node* getNode1(){return node1;}
};


typedef struct {
    double r;       // percent
    double g;       // percent
    double b;       // percent
} rgb;

typedef struct {
    double h;       // angle in degrees
    double s;       // percent
    double v;       // percent
} hsv;

static hsv   rgb2hsv(rgb in);
static rgb   hsv2rgb(hsv in);

hsv rgb2hsv(rgb in);
rgb hsv2rgb(hsv in);
//rgb:: negative: 230, 0, 0 to positive: 0, 0, 230,
//hsv:: negative: 0, 100%, 90% to positive: 240, 100%, 90%..
rgb getRGBValue(double startValue, double endValue, double value);
#endif /* graph_hpp */
