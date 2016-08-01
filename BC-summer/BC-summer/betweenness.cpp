//
//  betweenness.cpp
//  BC
//
//  Created by Anakin on 16/7/13.
//  Copyright © 2016年 Anakin. All rights reserved.
//

#include "betweenness.hpp"
using namespace boost;
using namespace boost::detail::graph;
using namespace std;
int Betweenness::compute(std::vector<Node *> &nodes, bool needDerivs)
{
    int n = (int)nodes.size();
    typedef adjacency_list<vecS, vecS, directedS, no_property, property<edge_weight_t, float, property<edge_weight2_t, float>>> Graph;
    Graph g;
    
    int numEdges = 0;
    vector<float> weights;
    vector<weighted_edge> tempEdges;
    
    for(int i = 0; i < n; i++)
    {
        for(unsigned int k = 0; k < nodes[i]->edges.size(); k++)
        {
            Edge* edge = nodes[i]->edges[k];
            add_edge(edge->getNode0()->getIndex(), edge->getNode1()->getIndex(), g);
            
            float weight = 1.0;
            //custom the weight value here!!!
            weights.push_back(weight);
            weighted_edge e;
            e.source = edge->getNode0()->getIndex();
            e.target = edge->getNode1()->getIndex();
            e.weight = weight;
            tempEdges.push_back(e);
            numEdges++;
        }
    }
    weighted_edge* edgesweighted = new weighted_edge[numEdges];
    for(int i = 0; i < numEdges; i++)
    {
        edgesweighted[i].source = tempEdges[i].source;
        edgesweighted[i].target = tempEdges[i].target;
        edgesweighted[i].weight = tempEdges[i].weight;
    }
    
    //the value of property_map, should be visited using iterator..
    property_map<Graph, edge_weight_t>::type w = get(edge_weight, g);
    float* weights0 = new float[numEdges];
    float* wp = weights0;
    for(int i = 0; i < numEdges; i++)
        wp[i] = weights[i];
    
    graph_traits<Graph>::edge_iterator e, e_end;
    for(boost::tie(e, e_end) = edges(g); e != e_end; ++e)
    {
        w[*e] = *wp++;
    }
    
    //create space for centralities..
    std::vector<double> edgeCentrality(num_edges(g));
    std::vector<double> vertex_centralities(num_vertices(g));
    std::vector<double> edge_centralities(num_edges(g));
    run_weighted_test((Digraph*)0, n, edgesweighted, numEdges, edge_centralities, vertex_centralities);
    
    
#ifdef DEBGU_CENTRALITY
    //output the results:
    cout << "edge centralities: " << endl;
    for(int i = 0; i < edge_centralities.size(); ++i)
    {
        cout << i << ":" <<edge_centralities[i] << " ";
    }
    cout << endl;
    cout << "vertex centralities: " << endl;
    for(int i = 0; i < vertex_centralities.size(); ++i)
    {
        cout << i << ":" << vertex_centralities[i] << " ";
    }
    cout << endl;
#endif
    //save the node centrality values:
    for(int i = 0; i < vertex_centralities.size(); ++i)
    {
        nodes[i]->centralityValue = vertex_centralities[i];
    }
//******************** for the sensitivity *****************************
    if(needDerivs)
    {
        for(int i = 0; i < n; i++)
        {
            int numEdges = 0;
            for(int gi = 0; gi < n; gi++)
            {
                float sumWeights = 0;
                for(unsigned int k = 0; k < nodes[gi]->edges.size(); k++)
                {
                    Edge* edge = nodes[gi]->edges[k];
                    float weight = 1.0;
                    //custom the weight value here!!!
                    sumWeights += weight;
                }
                for(unsigned int k=0;k<nodes[gi]->edges.size();k++) {
                    Edge * edge = nodes[gi]->edges[k];
                    int gj = edge->getNode1()->getIndex();
                    float weight = 1.0;
                    //custom the weight value here!!!
                    weights[numEdges] = ((i==gi || i==gj) && sumWeights!=0)? weight + weight/sumWeights: weight;
                    numEdges++;
                }
            }
            

            
            
            
            graph_traits < Graph >::edge_iterator e, e_end;
            int k=0;
            for (boost::tie(e, e_end) = edges(g); e != e_end; ++e) {
                //printf("Weight[%d] = %f\n", k, weights[k]);
                w[*e] = weights[k++];
            }
            std::vector<double> centrality2(n);
            
            brandes_betweenness_centrality(g,
                                           centrality_map(
                                                          make_iterator_property_map(centrality2.begin(),
                                                                                     get(vertex_index, g),
                                                                                     double()))
                                           .vertex_index_map(get(vertex_index, g))
                                           .weight_map(get(boost::edge_weight, g))
                                           );
            
#ifdef DEBUG_SENSITIVITY
            cout << "sensitivity: " << i << endl;
            for(int j = 0; j < centrality2.size(); ++j)
            {
                cout << j << ":" <<centrality2[j]-vertex_centralities[j] << " ";
            }
            cout << endl;
#endif
        
        
            for(int j = 0; j < centrality2.size(); ++j)
            {
                nodes[i]->sensitivityValues.push_back(centrality2[j]-vertex_centralities[j]);
            }
#ifdef DEBUG_CHECK_ITERATION
            std::cout<<"finish iteration: " << i << endl;
#endif
        }
        
        //compute the mean value and variance value
        std::map<int, double> sensitivityHash;
        std::map<int, double>::iterator iter;
        
        for(int i = 0; i < nodes.size(); i++)
            sensitivityHash.insert(make_pair(i, 0));
        
        for(int i = 0; i < nodes.size(); i++)
        {
            for(int j = 0; j < nodes[i]->sensitivityValues.size(); j++)
            {
                if(j != i)
                {
                    sensitivityHash[j] += nodes[i]->sensitivityValues[j];
                }
            }
        }
    
        for(int i = 0; i < nodes.size(); i++)
        {
            nodes[i]->sensitivityMean = sensitivityHash[i] / (nodes.size() - 1);
        }
        
        std::map<int, double> sensitivityHash2;
        
        for(int i = 0; i < nodes.size(); i++)
            sensitivityHash2.insert(make_pair(i, 0));
        
        for(int i = 0; i < nodes.size(); i++)
        {
            for(int j = 0; j < nodes[i]->sensitivityValues.size(); j++)
            {
                if(j != i)
                {
                    sensitivityHash2[j] += (nodes[i]->sensitivityValues[j]-nodes[j]->sensitivityMean)*(nodes[i]->sensitivityValues[j]-nodes[j]->sensitivityMean);
                }
            }
        }
        
        for(int i = 0; i < nodes.size(); i++)
            nodes[i]->sensitivityVariance = sqrt(sensitivityHash2[i]/(nodes.size()-1));
        
#ifdef DEBUG_MEAN_VARIANCE
        for(int i = 0; i < nodes.size(); i++)
        {
            cout << "the mean value and variance for node "<< i << endl;
            cout << nodes[i]->sensitivityMean << endl;
            cout << nodes[i]->sensitivityVariance << endl;
        }
#endif
#ifdef DEBUG_SENSITIVITY_CSV
        ofstream outputCSV;
        outputCSV.open("/Users/anakin/github/network-summer/BC-summer/BC-summer/sensitivityRes.csv");
        for(int j = 0; j < nodes.size(); ++j)
        {
            outputCSV << j << ",";
            for(int i = 0; i < nodes[j]->sensitivityValues.size(); i++)
            {
                outputCSV << nodes[j]->sensitivityValues[i] << ",";
            }
            outputCSV << "\n";
        }
        outputCSV.close();
#endif
        
    }

//****************************************************
    return 0;
}

template <typename GraphW>
void run_weighted_test(GraphW*, int V, weighted_edge edge_init[], int E, std::vector<double>& centralities, std::vector<double>& node_centralities)
{
    GraphW g(V);
    typedef typename graph_traits<GraphW>::vertex_descriptor Vertex;
    typedef typename graph_traits<GraphW>::vertex_iterator vertex_iterator;
    typedef typename graph_traits<GraphW>::edge_descriptor Edge;
    std::vector<Vertex> vertices(V);
    {
        vertex_iterator v, v_end;
        int index = 0;
        for(tie(v, v_end) = boost::vertices(g); v != v_end; ++v, ++index)
        {
            put(vertex_index, g, *v, index);
            vertices[index] = *v;
        }
    }
    std::vector<Edge> edges(E);
    for(int e = 0; e < E; ++e)
    {
        //return a pair, the first one is an edge descriptor, and the second one is a bool..
        edges[e] = add_edge(vertices[edge_init[e].source], vertices[edge_init[e].target], g).first;
        put(edge_weight, g, edges[e], edge_init[e].weight);
        put(edge_index, g, edges[e], e);
    }
    
    std::vector<double> centrality(V);
    std::vector<double> edge_centrality2(E);
    //***************** important *********************
    brandes_betweenness_centrality(g,
                                   centrality_map(make_iterator_property_map(centrality.begin(), get(vertex_index, g), double()))
                                   .edge_centrality_map(make_iterator_property_map(edge_centrality2.begin(), get(edge_index, g), double()))
                                   .vertex_index_map(get(vertex_index, g))
                                   .weight_map(get(edge_weight, g)));
    
    
    //***************** my try to read sigma, distance values from the function *****
    typedef typename graph_traits<GraphW>::degree_size_type degree_size_type;
    typedef typename graph_traits<GraphW>::edge_descriptor edge_descriptor;
    typedef decltype(make_iterator_property_map(centrality.begin(), get(vertex_index, g), double())) a_centrality_map;
    typedef typename property_traits<a_centrality_map>::value_type centrality_type;

    
    std::vector<std::vector<edge_descriptor> > incoming(V);
    std::vector<centrality_type> distance(V);
    std::vector<centrality_type> dependency(V);
    std::vector<degree_size_type> path_count(V);
    
//    auto incoming_map = make_iterator_property_map(incoming.begin(), get(vertex_index, g));
//    auto distance_map = make_iterator_property_map(distance.begin(), get(vertex_index, g));
//    auto dependency_map = make_iterator_property_map(dependency.begin(), get(vertex_index, g));
//    auto path_count_map = make_iterator_property_map(path_count.begin(), get(vertex_index, g));
    
    brandes_betweenness_centrality(g,
                                   make_iterator_property_map(centrality.begin(), get(vertex_index, g), double()),
                                   make_iterator_property_map(edge_centrality2.begin(), get(edge_index, g), double()),
                                   make_iterator_property_map(incoming.begin(), get(vertex_index, g)),
                                   make_iterator_property_map(distance.begin(), get(vertex_index, g)),
                                   make_iterator_property_map(dependency.begin(), get(vertex_index, g)),
                                   make_iterator_property_map(path_count.begin(), get(vertex_index, g)),
                                   get(vertex_index, g),
                                   get(edge_weight, g));
    
    
    auto aa = path_count[0];
    
    /*
     brandes_betweenness_centrality(const Graph& g,
     CentralityMap centrality,     // C_B
     EdgeCentralityMap edge_centrality_map,
     IncomingMap incoming, // P
     DistanceMap distance,         // d
     DependencyMap dependency,     // delta
     PathCountMap path_count,      // sigma
     VertexIndexMap vertex_index,
     WeightMap weight_map
     */

//********************************
    
    centralities.clear();
    for(int i = 0; i < E; i++)
    {
        centralities.push_back(edge_centrality2[i]);
    }
    node_centralities.clear();
    for(int i = 0; i < V; i++)
    {
        node_centralities.push_back(centrality[i]);
    }
}


void Betweenness::brandes_implementation(std::vector<Node*> &nodes)
{
    //global store..
    for(int i = 0; i < nodes.size(); i++)
    {
        std::map<Node*, double> cost_per_node;

        for(int j = 0; j < nodes[i]->edges.size(); j++)
        {
            cost_per_node[nodes[i]->edges[j]->node1] = nodes[i]->edges[j]->weight;
        }
        for(int j = 0; j < nodes.size(); j++)
        {
            if(cost_per_node.find(nodes[j]) != cost_per_node.end())
                continue;
            cost_per_node[nodes[j]] = DBL_MAX;
        }
        cost_per_node[nodes[i]] = 0;
        cost_store[nodes[i]] = cost_per_node;
//        for(auto iter = cost_per_node.begin(); iter != cost_per_node.end(); iter++)
//            cost_store[nodes[i]][iter->first] = iter->second;
    }

    
    for(int i = 0; i < nodes.size(); i++)
        CB[nodes[i]] = 0;
    for(int i = 0; i < nodes.size(); i++)
    {
        auto r = nodes[i];
        std::stack<Node*> S;
        std::queue<Node*> Q;
        std::map<Node*, std::vector<Node*> > P;
        std::map<Node*, int> sigma;
        std::map<Node*, double> distance;
        for(int j = 0; j < nodes.size(); j++)
        {
            sigma[nodes[j]] = 0;
            distance[nodes[j]] = DBL_MAX;
            if(j == i)
            {
                sigma[nodes[j]] = 1;
                distance[nodes[j]] = 0;
            }
        }
        
        Q.push(r);
        
        //BFS traversal
        while(!Q.empty())
        {
            auto v = Q.front();
            Q.pop();
            S.push(v);
            for(int e = 0; e < v->edges.size(); e++)
            {
                auto w = v->edges[e]->getNode1();
                if(distance[w] == DBL_MAX)
                {
                    Q.push(w);
                    distance[w] = distance[v] + 1;
                }
                if(distance[w] == distance[v] + 1)
                {
                    sigma[w] = sigma[w] + sigma[v];
                    P[w].push_back(v);
                }
            }
        }
        
        // dependency accumulation
        std::map<Node*, double> delta;
        for(int j = 0; j < nodes.size(); j++)
            delta[nodes[j]] = 0;
        while (!S.empty()) {
            auto w = S.top();
            S.pop();
            for(int j = 0; j < P[w].size(); j++)
            {
                auto v = P[w][j];
                delta[v] = delta[v] + (sigma[v]*1.0)/(sigma[w])*(1 + delta[w]);
            }
            if(w != r)
            {
                CB[w] = CB[w] + delta[w];
            }
        }
        
        P_store[r] = P;
        sigma_store[r] = sigma;
        distance_store[r] = distance;
    }
//    std::cout << "distance store " << distance_store[nodes[50]][nodes[43]] << std::endl;
//    std::cout << "cost store " << cost_store[nodes[50]][nodes[43]] << std::endl;

}
void Betweenness::insertEdge(Node* src, Node* dest, double cost)
{
    sigma_old.clear();
    distance_old.clear();
    trackLost.clear();
    pairsDone.clear();
    
    //mind this, undirected graph...
    cost_store[src][dest] = cost;
//    cost_store[dest][src] = cost;
//    std::cout << "cost store " << cost_store[src][dest] << std::endl;
    std::vector<Node*> sinks = insertUpdate(dest, src, src);
    std::vector<Node*> sources = insertUpdate(src, dest, dest);
    for(int i = 0; i < sinks.size(); i++)
        insertUpdate(src, dest, sinks[i]);
    for(int i = 0; i < sources.size(); i++)
        insertUpdate(dest, src, sources[i]);
    increaseBetweenness();
}
std::vector<Node*> Betweenness::insertUpdate(Node* src, Node* dest, Node* z)
{
    std::vector<std::pair<Node*, Node*> > workSet;
    std::vector<Node*> visitedVertices;
    std::vector<Node*> affectedVertices;
    
    workSet.push_back(make_pair(src, dest));
    visitedVertices.push_back(src);
    
    while (workSet.size() != 0) {
        Node* x = workSet.back().first;
        Node* y = workSet.back().second;
        workSet.pop_back();
        double alt = cost_store[x][y] + distance_store[y][z];
        if(alt < distance_store[x][z])
        {
            if(!isIn(make_pair(x, z), sigma_old))
            {
                distance_old[x][z] = distance_store[x][z];
                sigma_old[x][z] = sigma_store[x][z];
                reduceBetweenness(x, z);
                sigma_store[x][z] = 0;
                P_store[x][z].clear();
            }
            if(isIn(make_pair(x, z), pairsDone))
                pairsDone.erase(find(pairsDone.begin(), pairsDone.end(), make_pair(x, z)));
            distance_store[x][z] = alt;
        }
        if(alt == distance_store[x][z] && distance_store[x][z] != DBL_MAX)
        {
            if(!isIn(make_pair(x, z), pairsDone))
            {
                if(!isIn(make_pair(x, z), sigma_old))
                    reduceBetweenness(x, z);
                if(sigma_store[x][z] != 0)
                    sigma_old[x][z] = sigma_store[x][z];
                
                sigma_store[x][z] = sigma_store[x][z] + (sigma_store[x][src] * 1 * sigma_store[dest][z]);
                P_store[x][y].push_back(x);
                for(int iter = 0; iter < P_store[y][z].size(); iter++)
                {
                    P_store[x][z].push_back(P_store[y][z][iter]);
                }
            }
            pairsDone.push_back(make_pair(x, z));
            affectedVertices.push_back(x);
            // need pred(x)...
            for(auto iter = x->pred.begin(); iter != x->pred.end(); iter++)
            {
                Node* u = iter->second;
                if(SP(u, x, src) && std::find(visitedVertices.begin(), visitedVertices.end(), u) != visitedVertices.end())
                {
                    workSet.push_back(make_pair(u, x));
                    visitedVertices.push_back(u);
                }
            }
        }
    }
    return affectedVertices;
}
void Betweenness::reduceBetweenness(Node* a, Node* z)
{
    if(sigma_old[a][z] == 0)
        return;
    std::vector<Node*> known;
    std::vector<Node*> stack;
    for(int i = 0; i < P_store[a][z].size(); i++)
    {
        Node* n = P_store[a][z][i];
        if(distance_store[a][z] != distance_old[a][n] + distance_old[n][z])
            continue;
        else if(a != n && n != z)
        {
            CB[n] = CB[n] - (sigma_old[a][n] * sigma_old[n][z] / sigma_old[a][z]);
            trackLost.push_back(make_tuple(a, z, n));
        }
        stack.push_back(n);
        known.push_back(n);
    }
    while (stack.size() != 0) {
        Node* p = stack.back();
        stack.pop_back();
        known.push_back(p);
        for(int i = 0; i < P_store[a][p].size(); i++)
        {
            Node* n = P_store[a][p][i];
            if(distance_store[a][z] != distance_old[a][n] + distance_old[n][z])
                continue;
            else if(a != n && n != z && std::find(known.begin(), known.end(), n) != known.end())
            {
                CB[n] = CB[n] - (sigma_old[a][n] * sigma_old[n][z] / sigma_old[a][z]);
                trackLost.push_back(make_tuple(a, z, n));
            }
            stack.push_back(n);
            known.push_back(n);
        }
    }
    std::vector<Node*> alreadyKnown;
    for(int i = 0; i < known.size(); i++)
        alreadyKnown.push_back(known[i]);
    alreadyKnown.push_back(a);
    for(int i = 0; i < trackLost.size(); i++)
    {
        Node* v = std::get<2>(trackLost[i]);
        Node* v1 = std::get<0>(trackLost[i]);
        Node* v2 = std::get<1>(trackLost[i]);
        if(std::find(known.begin(), known.end(), v1) != known.end() && std::find(known.begin(), known.end(), v2) != known.end()
           && distance_store[a][z] == distance_old[a][v] + distance_old[v][z])
        {
            if(std::find(alreadyKnown.begin(), alreadyKnown.end(), v) != alreadyKnown.end())
            {
                CB[v] = CB[v] - (sigma_old[a][v] * sigma_old[v][z] / sigma_old[a][z]);
                alreadyKnown.push_back(v);
                trackLost.push_back(std::make_tuple(a, z, v));
            }
        }
    }
    
}
void Betweenness::increaseBetweenness()
{
    for(auto outer = sigma_old.begin(); outer != sigma_old.end(); outer++)
    {
        auto src = outer->first;
        for(auto inner = outer->second.begin(); inner != outer->second.end(); inner++)
        {
            auto dest = inner->first;
            std::vector<Node*> known;
            std::vector<Node*> stack;
            for(int i = 0; i < P_store[src][dest].size(); i++)
            {
                Node* n = P_store[src][dest][i];
                stack.push_back(n);
                known.push_back(n);
                if(src != n && n != dest)
                    CB[n] = CB[n] + (sigma_store[src][n] * sigma_store[n][dest] / sigma_store[src][dest]);
            }
            while (stack.size() != 0) {
                Node* n = stack.back();
                stack.pop_back();
                for(int i = 0; i < P_store[src][n].size(); i++)
                {
                    Node* p = P_store[src][n][i];
                    if(p != src && p != dest && std::find(known.begin(), known.end(), p) == known.end())
                    {
                        stack.push_back(p);
                        known.push_back(p);
                         CB[p] = CB[p] + (sigma_store[src][p] * sigma_store[p][dest] / sigma_store[src][dest]);
                    }
                }
            }
        }
    }
}

bool Betweenness::isIn(std::pair<Node*, Node*> key, std::map<Node*, std::map<Node*, double> > container)
{
    if(container.find(key.first) != container.end())
        if(container[key.first].find(key.second) != container[key.first].end())
            return true;
    return false;
}
bool Betweenness::isIn(std::pair<Node*, Node*> key, std::map<Node*, std::map<Node*, int> > container)
{
    if(container.find(key.first) != container.end())
        if(container[key.first].find(key.second) != container[key.first].end())
            return true;
    return false;
}
bool Betweenness::isIn(std::pair<Node*, Node*> value, std::vector<std::pair<Node*, Node*> > container)
{
    if(find(container.begin(), container.end(), value) != container.end())
        return true;
    return false;
}
bool Betweenness::isIn(std::tuple<Node*, Node*, Node*> value, std::vector<std::tuple<Node*, Node*, Node*> > container)
{
    if(find(container.begin(), container.end(), value) != container.end())
        return true;
    return false;
}
bool Betweenness::SP(Node *x, Node *y, Node *z)
{
    if(isIn(make_pair(x, z), distance_store))
    {
        if(distance_store[x][z] != DBL_MAX)
        {
            if(distance_store[x][z] == cost_store[x][y] + distance_store[y][z])
                return true;
        }
    }
    return false;
}









