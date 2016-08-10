//
//  betweenness.hpp
//  BC
//
//  Created by Anakin on 16/7/13.
//  Copyright © 2016年 Anakin. All rights reserved.
//

#ifndef betweenness_hpp
#define betweenness_hpp

#include <iostream>
#include <utility>                   // for std::pair
#include <algorithm>                 // for std::for_each
#include <limits.h>
#include <vector>
#include <map>
#include <queue>
#include <stack>
#include <set>
#include <list>
#include <tuple>
#include <fstream>
#include <float.h>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/betweenness_centrality.hpp>
#include <boost/graph/bc_clustering.hpp>
#include <boost/graph/overloading.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/relax.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/if.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/named_function_params.hpp>

#include "graph.h"
#include "define.h"
#include "timelogger.h"



using namespace boost;
struct weighted_edge
{
    int source, target;
    double weight;
};


typedef property<edge_weight_t, double, property<edge_index_t, std::size_t>> EdgeProperties;
typedef adjacency_list<listS, listS, undirectedS, property<vertex_index_t, int>, EdgeProperties> GraphW;
typedef adjacency_list<listS, listS, directedS, property<vertex_index_t, int>, EdgeProperties> Digraph;

template <typename GraphW>
void run_weighted_test(GraphW*, int V, weighted_edge edge_init[], int E, std::vector<double>& centralities, std::vector<double>& node_centralities);
template <typename GraphW>
void run_unweighted_test(GraphW*, int V, weighted_edge edge_init[], int E, std::vector<double>&centralities);



class Betweenness
{
public:
//    double btwValue[];
    std::unordered_map<Node*, double> CB;
    std::unordered_map<Node*, double> CB_transpose;

    std::unordered_map<Node*, std::unordered_map<Node*, double> >sensitivityRes;

    std::unordered_map<Node*, std::unordered_map<Node*, std::vector<Node*> > > P_store;
    std::unordered_map<Node*, std::unordered_map<Node*, int> > sigma_store;
    std::unordered_map<Node*, std::unordered_map<Node*, double> > distance_store;
    std::unordered_map<Node*, std::unordered_map<Node*, double> > cost_store;

    std::unordered_map<Node*, std::unordered_map<Node*, std::vector<Node*> > > P_store_transpose;
    std::unordered_map<Node*, std::unordered_map<Node*, int> > sigma_store_transpose;
    std::unordered_map<Node*, std::unordered_map<Node*, double> > distance_store_transpose;
    std::unordered_map<Node*, std::unordered_map<Node*, double> > cost_store_transpose;

    std::unordered_map<Node*, std::unordered_map<Node*, int> > sigma_old;
    std::unordered_map<Node*, std::unordered_map<Node*, double> > distance_old;
    std::vector<std::pair<Node*, Node*> > pairsDone;
    std::vector<std::tuple<Node*, Node*, Node*> > trackLost;




    int compute(std::vector<Node*> &nodes, bool needDerivs);
    void brandes_implementation_init(std::vector<Node*> &nodes);
    void brandes_implementation(std::vector<Node*> &nodes);
    void brandes_implementation_weighted(std::vector<Node*> &nodes, bool isTranspose);
    //void brandes_implementation_weighted_transpose(std::vector<Node*> &nodes);

    void insertEdge(Node* src, Node* dest, double cost);
    std::vector<Node*> insertUpdate(Node* src, Node* dest, Node* z, bool isTranspose);
    void reduceBetweenness(Node* a, Node* z, bool isTranspose);
    void increaseBetweenness(bool isTranspose);

    void calSensitivity(std::vector<Node*> &nodes, std::unordered_map<Node*, double> vertex_centralities);

//    void deleteEdge(Node* src, Node* dest);
//    std::vector<Node*> deleteUpdate(Node* src, Node* dest, Node* z);


    bool isIn(std::pair<Node*, Node*> key, std::unordered_map<Node*, std::unordered_map<Node*, double> > &container);
    bool isIn(std::pair<Node*, Node*> key, std::unordered_map<Node*, std::unordered_map<Node*, int> > &container);
    bool isIn(std::pair<Node*, Node*> value, std::vector<std::pair<Node*, Node*> > &container);
    bool isIn(std::tuple<Node*, Node*, Node*> value, std::vector<std::tuple<Node*, Node*, Node*> > &container);

    double getDistanceOldVal(Node* x, Node* y, bool isTranspose);
    int getSigmaOldVal(Node* x, Node* y, bool isTranspose);

    bool SP(Node* x, Node* y, Node* z, bool isTranspose);
};
#endif /* betweenness_hpp */
