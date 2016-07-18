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

#include "graph.hpp"
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
    int compute(std::vector<Node*> &nodes, bool needDerivs);
};
#endif /* betweenness_hpp */
