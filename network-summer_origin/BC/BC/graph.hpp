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
//#include "properties.hpp"
class Edge;
class Node
{
public:
    std::vector<Edge*> edges;
    int index;
//    std::string property;
    Node(){ index = -1; }
    ~Node()
    {
        for(int i = 0; i < edges.size(); i++)
            delete edges[i];
    }
    
    Edge* addEdge(Node* n);
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
    Node* node0;
    Node* node1;
    Edge(){ node0 = 0; node1 = 0; multiplicity = 0; weight = 0; }
    void setNode0(Node* n){ node0 = n; }
    void setNode1(Node* n){ node1 = n; }
    Node* getNode0(){return node0;}
    Node* getNode1(){return node1;}
};

#endif /* graph_hpp */
