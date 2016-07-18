//
//  graph.cpp
//  BC
//
//  Created by Anakin on 16/7/14.
//  Copyright © 2016年 Anakin. All rights reserved.
//

#include "graph.hpp"
#include <math.h>
#include <algorithm>
bool Node::operator==(Node &b)
{
    return index == b.index;
}
bool Node::operator==(int idx)
{
    return index == idx;
}

int Node::getIndex()
{
    return index;
}
void Node::setIndex(int idx)
{
    index = idx;
}

int Node::findEdge(Node *n)
{
    for(int i = 0; i < edges.size(); i++)
    {
        if(edges[i]->getNode1() == n)
        {
            return i;
        }
    }
    return -1;
}

Edge* Node::getEdge(int idx)
{
    for(int i = 0; i < edges.size(); i++)
    {
        if(edges[i]->getNode1()->getIndex() == idx)
        {
            return edges[i];
        }
    }
    return 0;
}
Edge* Node::addEdge(Node *n)
{
    int k = findEdge(n);
    if(k == -1)
    {
        Edge* e = new Edge;
        e->setNode0(this);
        e->setNode1(n);
        e->multiplicity = 1;
        edges.push_back(e);
        return e;
    }
    else
    {
        edges[k]->multiplicity = edges[k]->multiplicity + 1;
        return edges[k];
    }
}













