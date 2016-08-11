//
//  graph.cpp
//  BC
//
//  Created by Anakin on 16/7/14.
//  Copyright © 2016年 Anakin. All rights reserved.
//

#include "graph.h"
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
Edge* Node::addEdge(Node *n, double weight)
{
    int k = findEdge(n);
    if(k == -1)
    {
        Edge* e = new Edge;
        e->setNode0(this);
        e->setNode1(n);
        e->multiplicity = 1;
        e->weight = weight;      //needs to be modified..
        edges.push_back(e);
        n->pred[weight].push_back(this);

        Edge* e_transpose = new Edge;
        e_transpose->setNode0(n);
        e_transpose->setNode1(this);
        e_transpose->multiplicity = 1;
        e_transpose->weight = weight;      //needs to be modified..
        n->edges_transpose.push_back(e_transpose);
        pred_transpose[weight].push_back(n);


        return e;
    }
    else
    {
        edges[k]->multiplicity = edges[k]->multiplicity + 1;

        n->edges_transpose[k]->multiplicity = n->edges_transpose[k]->multiplicity + 1;

        return edges[k];
    }

}

hsv rgb2hsv(rgb in)
{
    hsv         out;
    double      min, max, delta;

    min = in.r < in.g ? in.r : in.g;
    min = min  < in.b ? min  : in.b;

    max = in.r > in.g ? in.r : in.g;
    max = max  > in.b ? max  : in.b;

    out.v = max;                                // v
    delta = max - min;
    if (delta < 0.00001)
    {
        out.s = 0;
        out.h = 0; // undefined, maybe nan?
        return out;
    }
    if( max > 0.0 ) { // NOTE: if Max is == 0, this divide would cause a crash
        out.s = (delta / max);                  // s
    } else {
        // if max is 0, then r = g = b = 0
        // s = 0, v is undefined
        out.s = 0.0;
        out.h = NAN;                            // its now undefined
        return out;
    }
    if( in.r >= max )                           // > is bogus, just keeps compilor happy
        out.h = ( in.g - in.b ) / delta;        // between yellow & magenta
    else
        if( in.g >= max )
            out.h = 2.0 + ( in.b - in.r ) / delta;  // between cyan & yellow
        else
            out.h = 4.0 + ( in.r - in.g ) / delta;  // between magenta & cyan

    out.h *= 60.0;                              // degrees

    if( out.h < 0.0 )
        out.h += 360.0;

    return out;
}

rgb hsv2rgb(hsv in)
{
    double      hh, p, q, t, ff;
    long        i;
    rgb         out;

    if(in.s <= 0.0) {       // < is bogus, just shuts up warnings
        out.r = in.v;
        out.g = in.v;
        out.b = in.v;
        return out;
    }
    hh = in.h;
    if(hh >= 360.0) hh = 0.0;
    hh /= 60.0;
    i = (long)hh;
    ff = hh - i;
    p = in.v * (1.0 - in.s);
    q = in.v * (1.0 - (in.s * ff));
    t = in.v * (1.0 - (in.s * (1.0 - ff)));

    switch(i) {
        case 0:
            out.r = in.v;
            out.g = t;
            out.b = p;
            break;
        case 1:
            out.r = q;
            out.g = in.v;
            out.b = p;
            break;
        case 2:
            out.r = p;
            out.g = in.v;
            out.b = t;
            break;

        case 3:
            out.r = p;
            out.g = q;
            out.b = in.v;
            break;
        case 4:
            out.r = t;
            out.g = p;
            out.b = in.v;
            break;
        case 5:
        default:
            out.r = in.v;
            out.g = p;
            out.b = q;
            break;
    }
    return out;
}

rgb getRGBValue(double startValue, double endValue, double value)
{
    rgb rgbNegativeMax;
    rgbNegativeMax.r = 230;
    rgbNegativeMax.g = 0;
    rgbNegativeMax.b = 0;

    rgb rgbPositiveMax;
    rgbPositiveMax.r = 0;
    rgbPositiveMax.g = 0;
    rgbPositiveMax.b = 230;

    hsv hsvNegativeMax = rgb2hsv(rgbNegativeMax);
    hsv hsvPositiveMax = rgb2hsv(rgbPositiveMax);

    double scale;
    rgb result;

    if(startValue < 0)
    {
        if(value < 0)
        {
            scale = log10(-value+1)/log10(-startValue+1);
            hsvNegativeMax.s *= scale;
            return hsv2rgb(hsvNegativeMax);
        }
        else
        {
            scale = log10(value+1)/log10(endValue+1);
            hsvPositiveMax.s *= scale;
            return hsv2rgb(hsvPositiveMax);
        }
    }
    else
    {
        scale = log10(value+1)/log10(endValue+1);
        hsvPositiveMax.s *= scale;
        return hsv2rgb(hsvPositiveMax);
    }
}
