//
//  betweenness.cpp
//  BC
//
//  Created by Anakin on 16/7/13.
//  Copyright © 2016年 Anakin. All rights reserved.
//

#include "betweenness.h"
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
                    //custom the weight value here!!!
                    sumWeights += edge->weight;
                }
                for(unsigned int k=0;k<nodes[gi]->edges.size();k++) {
                    Edge * edge = nodes[gi]->edges[k];
                    int gj = edge->getNode1()->getIndex();
                    float weight = edge->weight;
                    //custom the weight value here!!!
//                    cout << gi << " " << gj << " " << numEdges <<endl;
//                    cout << (((i==gi || i==gj) && sumWeights!=0)? weight - weight/sumWeights: weight) << endl;

                    weights[numEdges] = ((i==gi || i==gj) && sumWeights!=0)? weight - weight/sumWeights: weight;
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
                cout << j << ":" <<vertex_centralities[j] - centrality2[j] << " ";
            }
            cout << endl;
#endif


            for(int j = 0; j < centrality2.size(); ++j)
            {
                nodes[i]->sensitivityValues.push_back(vertex_centralities[j] - centrality2[j]);
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
void Betweenness::brandes_implementation_init(std::vector<Node*> &nodes)
{
    //global store..
    for(int i = 0; i < nodes.size(); i++)
    {
        std::unordered_map<Node*, double> cost_per_node;

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
    {
        std::unordered_map<Node*, double> cost_per_node;

        for(int j = 0; j < nodes[i]->edges_transpose.size(); j++)
        {
            cost_per_node[nodes[i]->edges_transpose[j]->node1] = nodes[i]->edges_transpose[j]->weight;
        }
        for(int j = 0; j < nodes.size(); j++)
        {
            if(cost_per_node.find(nodes[j]) != cost_per_node.end())
                continue;
            cost_per_node[nodes[j]] = DBL_MAX;
        }
        cost_per_node[nodes[i]] = 0;
        cost_store_transpose[nodes[i]] = cost_per_node;
    }
}

void Betweenness::brandes_implementation(std::vector<Node*> &nodes)
{
    //global store..
//    for(int i = 0; i < nodes.size(); i++)
//    {
//        std::unordered_map<Node*, double> cost_per_node;

//        for(int j = 0; j < nodes[i]->edges.size(); j++)
//        {
//            cost_per_node[nodes[i]->edges[j]->node1] = nodes[i]->edges[j]->weight;
//        }
//        for(int j = 0; j < nodes.size(); j++)
//        {
//            if(cost_per_node.find(nodes[j]) != cost_per_node.end())
//                continue;
//            cost_per_node[nodes[j]] = DBL_MAX;
//        }
//        cost_per_node[nodes[i]] = 0;
//        cost_store[nodes[i]] = cost_per_node;

//    }


    for(int i = 0; i < nodes.size(); i++)
        CB[nodes[i]] = 0;
    for(int i = 0; i < nodes.size(); i++)
    {
        auto r = nodes[i];
        std::stack<Node*> S;
        std::queue<Node*> Q;
        std::unordered_map<Node*, std::vector<Node*> > P;
        std::unordered_map<Node*, int> sigma;
        std::unordered_map<Node*, double> distance;
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
struct comparator {
    bool operator()(std::tuple<double, int, Node *, Node *> i, std::tuple<double, int, Node *, Node *> j) {
        if(std::get<0>(i) == std::get<0>(j))
            return std::get<1>(i) > std::get<1>(j);
        return std::get<0>(i) > std::get<0>(j);
    }
};


void Betweenness::brandes_implementation_weighted(std::vector<Node*> &nodes, bool isTranspose)
{

    std::unordered_map<Node*, std::unordered_map<Node*, std::vector<Node*> > > &P_now = isTranspose? P_store_transpose:P_store;
    std::unordered_map<Node*, std::unordered_map<Node*, int> > &sigma_now = isTranspose? sigma_store_transpose:sigma_store;
    std::unordered_map<Node*, std::unordered_map<Node*, double> > &distance_now = isTranspose? distance_store_transpose:distance_store;
    std::unordered_map<Node*, std::unordered_map<Node*, double> > &cost_now = isTranspose? cost_store_transpose:cost_store;
    std::unordered_map<Node*, double> &CB_now = isTranspose? CB_transpose:CB;

    for(int i = 0; i < nodes.size(); i++)
        CB_now[nodes[i]] = 0.0;
    for(int i = 0; i < nodes.size(); i++)
    {
        auto r = nodes[i];
        std::stack<Node*> S;
        std::unordered_map<Node*, std::vector<Node*> > P;
        std::unordered_map<Node*, int> sigma;
        std::unordered_map<Node*, double> distance;
        for(int j = 0; j < nodes.size(); j++)
        {
            sigma[nodes[j]] = 0;
            distance[nodes[j]] = DBL_MAX;
            if(j == i)
            {
                sigma[nodes[j]] = 1;
                //                distance[nodes[j]] = 0;
            }
        }

        std::unordered_map<Node*, double> seen;
        seen[r] = 0.0;
        int c = 0;
        priority_queue<std::tuple<double, int, Node *, Node *>, std::vector<std::tuple<double, int, Node *, Node *>>, comparator> Q;
//        std::queue<std::tuple<double, int, Node *, Node *>> Q;

        Q.push(make_tuple(0.0, c++, r, r));


        //dijkstra
        while(!Q.empty())
        {
            double dist = std::get<0>(Q.top());
            //            auto _ =std::get<1>(Q.front());
            auto pred = std::get<2>(Q.top());
            auto v = std::get<3>(Q.top());
            Q.pop();

            if(distance[v] != DBL_MAX)
                continue;
            sigma[v] += sigma[pred];

//            cout << "v: " << v->getIndex() << endl;
//            for(auto aa = sigma.begin(); aa != sigma.end(); aa++)
//            {
//                cout << aa->second << " ";
//            }
//            cout << endl;

            S.push(v);
            distance[v] = dist;

            int edgeSize = isTranspose? v->edges_transpose.size():v->edges.size();
            for(int e = 0; e < edgeSize; e++)
            {
                Edge* edge_now = isTranspose? v->edges_transpose[e] : v->edges[e];
                auto w = edge_now->getNode1();
                double vw_dist = dist + cost_now[edge_now->getNode0()][edge_now->getNode1()];//edge_now->weight;
//                cout << "v: " << v->getIndex() << " w: " << w->getIndex() << endl;
//                for(auto aa = sigma.begin(); aa != sigma.end(); aa++)
//                {
//                    cout << aa->second << " ";
//                }
//                cout << endl;


                if(distance[w] == DBL_MAX && (seen.find(w) == seen.end() || vw_dist < seen[w]))
                {
                    seen[w] = vw_dist;
                    Q.push(make_tuple(vw_dist, c++, v, w));
                    sigma[w] = 0;
                    P[w].clear();
                    P[w].push_back(v);
                }
                else if(vw_dist == seen[w])
                {
                    sigma[w] = sigma[w] + sigma[v];
                    P[w].push_back(v);
                }
            }
        }


        //why this happens??
        for(auto sigma_iter = sigma.begin(); sigma_iter != sigma.end(); sigma_iter++)
            sigma_iter->second /= 2;



        // dependency accumulation
        std::map<Node*, double> delta;
        for(int j = 0; j < nodes.size(); j++)
            delta[nodes[j]] = 0;
        while (!S.empty()) {
            auto w = S.top();
            S.pop();
            double coeff = (1.0 + delta[w]) / sigma[w];
            for(int j = 0; j < P[w].size(); j++)
            {
                auto v = P[w][j];
//                delta[v] = delta[v] + (double)(sigma[v]*1.0)/(sigma[w])*(1 + delta[w]);
                delta[v] += sigma[v] * coeff;
                //cout << "v:" << v->getIndex() << " delta[v]: " << delta[v] << endl;
            }
            if(w != r)
            {
                CB_now[w] = CB_now[w] + delta[w];
                //cout << "w:" << w->getIndex() << " CB[w]: " << CB[w] << endl;
            }
        }

        P_now[r] = P;
        sigma_now[r] = sigma;
        distance_now[r] = distance;
    }
}

//void Betweenness::brandes_implementation_weighted_transpose(std::vector<Node*> &nodes)
//{

//    for(int i = 0; i < nodes.size(); i++)
//        CB_transpose[nodes[i]] = 0.0;
//    for(int i = 0; i < nodes.size(); i++)
//    {
//        auto r = nodes[i];
//        std::stack<Node*> S;
//        std::unordered_map<Node*, std::vector<Node*> > P;
//        std::unordered_map<Node*, int> sigma;
//        std::unordered_map<Node*, double> distance;
//        for(int j = 0; j < nodes.size(); j++)
//        {
//            sigma[nodes[j]] = 0;
//            distance[nodes[j]] = DBL_MAX;
//            if(j == i)
//            {
//                sigma[nodes[j]] = 1;
//                //                distance[nodes[j]] = 0;
//            }
//        }

//        std::unordered_map<Node*, double> seen;
//        seen[r] = 0.0;
//        int c = 0;
//        priority_queue<std::tuple<double, int, Node *, Node *>, std::vector<std::tuple<double, int, Node *, Node *>>, comparator> Q;
////        std::queue<std::tuple<double, int, Node *, Node *>> Q;

//        Q.push(make_tuple(0.0, c++, r, r));


//        //dijkstra
//        while(!Q.empty())
//        {
//            double dist = std::get<0>(Q.top());
//            //            auto _ =std::get<1>(Q.front());
//            auto pred = std::get<2>(Q.top());
//            auto v = std::get<3>(Q.top());
//            Q.pop();

//            if(distance[v] != DBL_MAX)
//                continue;
//            sigma[v] += sigma[pred];

////            cout << "v: " << v->getIndex() << endl;
////            for(auto aa = sigma.begin(); aa != sigma.end(); aa++)
////            {
////                cout << aa->second << " ";
////            }
////            cout << endl;

//            S.push(v);
//            distance[v] = dist;

//            for(int e = 0; e < v->edges.size(); e++)
//            {
//                auto w = v->edges[e]->getNode1();
//                double vw_dist = dist + v->edges[e]->weight;
////                cout << "v: " << v->getIndex() << " w: " << w->getIndex() << endl;
////                for(auto aa = sigma.begin(); aa != sigma.end(); aa++)
////                {
////                    cout << aa->second << " ";
////                }
////                cout << endl;


//                if(distance[w] == DBL_MAX && (seen.find(w) == seen.end() || vw_dist < seen[w]))
//                {
//                    seen[w] = vw_dist;
//                    Q.push(make_tuple(vw_dist, c++, v, w));
//                    sigma[w] = 0;
//                    P[w].clear();
//                    P[w].push_back(v);
//                }
//                else if(vw_dist == seen[w])
//                {
//                    sigma[w] = sigma[w] + sigma[v];
//                    P[w].push_back(v);
//                }
//            }
//        }


//        //why this happens??
//        for(auto sigma_iter = sigma.begin(); sigma_iter != sigma.end(); sigma_iter++)
//            sigma_iter->second /= 2;



//        // dependency accumulation
//        std::map<Node*, double> delta;
//        for(int j = 0; j < nodes.size(); j++)
//            delta[nodes[j]] = 0;
//        while (!S.empty()) {
//            auto w = S.top();
//            S.pop();
//            double coeff = (1.0 + delta[w]) / sigma[w];
//            for(int j = 0; j < P[w].size(); j++)
//            {
//                auto v = P[w][j];
////                delta[v] = delta[v] + (double)(sigma[v]*1.0)/(sigma[w])*(1 + delta[w]);
//                delta[v] += sigma[v] * coeff;
//                //cout << "v:" << v->getIndex() << " delta[v]: " << delta[v] << endl;
//            }
//            if(w != r)
//            {
//                CB_transpose[w] = CB_transpose[w] + delta[w];
//                //cout << "w:" << w->getIndex() << " CB[w]: " << CB[w] << endl;
//            }
//        }

//        P_store_transpose[r] = P;
//        sigma_store_transpose[r] = sigma;
//        distance_store_transpose[r] = distance;
//    }
//}

void Betweenness::calSensitivityIncremental(std::vector<Node*> &nodes, std::unordered_map<Node*, double> vertex_centralities)
{

    //***********************
    //need to store:
    /*
        std::unordered_map<Node*, std::unordered_map<Node*, std::vector<Node*> > > P_store;
        std::unordered_map<Node*, std::unordered_map<Node*, int> > sigma_store;
        std::unordered_map<Node*, std::unordered_map<Node*, double> > distance_store;
        std::unordered_map<Node*, std::unordered_map<Node*, double> > cost_store;

        std::unordered_map<Node*, std::unordered_map<Node*, std::vector<Node*> > > P_store_transpose;
        std::unordered_map<Node*, std::unordered_map<Node*, int> > sigma_store_transpose;
        std::unordered_map<Node*, std::unordered_map<Node*, double> > distance_store_transpose;
        std::unordered_map<Node*, std::unordered_map<Node*, double> > cost_store_transpose;
     */
    std::unordered_map<Node*, std::unordered_map<Node*, std::vector<Node*> > > P_store_;
    std::unordered_map<Node*, std::unordered_map<Node*, int> > sigma_store_;
    std::unordered_map<Node*, std::unordered_map<Node*, double> > distance_store_;
    std::unordered_map<Node*, std::unordered_map<Node*, double> > cost_store_;

    std::unordered_map<Node*, std::unordered_map<Node*, std::vector<Node*> > > P_store_transpose_;
    std::unordered_map<Node*, std::unordered_map<Node*, int> > sigma_store_transpose_;
    std::unordered_map<Node*, std::unordered_map<Node*, double> > distance_store_transpose_;
    std::unordered_map<Node*, std::unordered_map<Node*, double> > cost_store_transpose_;

    for(auto iter = P_store.begin(); iter != P_store.end(); iter++)
    {
        std::unordered_map<Node*, std::vector<Node*> > P;
        for(auto iter_inner = iter->second.begin(); iter_inner != iter->second.end(); iter_inner++)
        {
            for(int ii = 0; ii < iter_inner->second.size(); ii++)
                P[iter_inner->first].push_back(iter_inner->second[ii]);
        }
        P_store_[iter->first] = P;
    }
    for(auto iter = P_store_transpose.begin(); iter != P_store_transpose.end(); iter++)
    {
        std::unordered_map<Node*, std::vector<Node*> > P;
        for(auto iter_inner = iter->second.begin(); iter_inner != iter->second.end(); iter_inner++)
        {
            for(int ii = 0; ii < iter_inner->second.size(); ii++)
                P[iter_inner->first].push_back(iter_inner->second[ii]);
        }
        P_store_transpose_[iter->first] = P;
    }
    for(auto iter = sigma_store.begin(); iter != sigma_store.end(); iter++)
    {
        std::unordered_map<Node*, int > sigma;
        for(auto iter_inner = iter->second.begin(); iter_inner != iter->second.end(); iter_inner++)
        {
                sigma[iter_inner->first] = iter_inner->second;
        }
        sigma_store_[iter->first] = sigma;
    }
    for(auto iter = sigma_store_transpose.begin(); iter != sigma_store_transpose.end(); iter++)
    {
        std::unordered_map<Node*, int > sigma;
        for(auto iter_inner = iter->second.begin(); iter_inner != iter->second.end(); iter_inner++)
        {
                sigma[iter_inner->first] = iter_inner->second;
        }
        sigma_store_transpose_[iter->first] = sigma;
    }
    for(auto iter = distance_store.begin(); iter != distance_store.end(); iter++)
    {
        std::unordered_map<Node*, double > distance;
        for(auto iter_inner = iter->second.begin(); iter_inner != iter->second.end(); iter_inner++)
        {
                distance[iter_inner->first] = iter_inner->second;
        }
        distance_store_[iter->first] = distance;
    }
    for(auto iter = distance_store_transpose.begin(); iter != distance_store_transpose.end(); iter++)
    {
        std::unordered_map<Node*, double > distance;
        for(auto iter_inner = iter->second.begin(); iter_inner != iter->second.end(); iter_inner++)
        {
                distance[iter_inner->first] = iter_inner->second;
        }
        distance_store_transpose_[iter->first] = distance;
    }
    for(auto iter = cost_store.begin(); iter != cost_store.end(); iter++)
    {
        std::unordered_map<Node*, double > cost;
        for(auto iter_inner = iter->second.begin(); iter_inner != iter->second.end(); iter_inner++)
        {
                cost[iter_inner->first] = iter_inner->second;
        }
        cost_store_[iter->first] = cost;
    }
    for(auto iter = cost_store_transpose.begin(); iter != cost_store_transpose.end(); iter++)
    {
        std::unordered_map<Node*, double > cost;
        for(auto iter_inner = iter->second.begin(); iter_inner != iter->second.end(); iter_inner++)
        {
                cost[iter_inner->first] = iter_inner->second;
        }
        cost_store_transpose_[iter->first] = cost;
    }



    //***********************
        int n = nodes.size();
        for(int i = 0; i < n; i++)
        {

//            brandes_implementation_init(nodes);
//            brandes_implementation_weighted(nodes, false);
//            brandes_implementation_weighted(nodes, true);

            for(int gi = 0; gi < n; gi++)
            {
                float sumWeights = 0;
                for(unsigned int k = 0; k < nodes[gi]->edges.size(); k++)
                {
                    Edge* edge = nodes[gi]->edges[k];
//                    float weight = 1.0;
                    //custom the weight value here!!!
                    sumWeights += edge->weight;
                }
                for(unsigned int k=0;k<nodes[gi]->edges.size();k++) {
                    Edge * edge = nodes[gi]->edges[k];
                    int gj = edge->getNode1()->getIndex();
                    float weight = edge->weight;
                    //custom the weight value here!!!
                    if((i==gi || i==gj) && sumWeights!=0)
                    {
//                        edge->weight = weight - weight/sumWeights;
//                        cout << edge->getNode0()->getIndex() << " " << edge->getNode1()->getIndex() << endl;
//                        cout << weight - weight/sumWeights << endl;
                        insertEdge(edge->getNode0(), edge->getNode1(), weight - weight/sumWeights);
                    }
                    //edge->weight = ((i==gi || i==gj) && sumWeights!=0)? weight + weight/sumWeights: weight;
                }
            }

#ifdef DEBUG_SENSITIVITY
            cout << "sensitivity incremental: " << i << endl;
//            for(auto j = CB.begin(); j != CB.end(); j++)
//            {
//                cout << j->first->getIndex() << ":" << vertex_centralities[j->first] - j->second << " ";
//            }
            for(int j = 0; j < nodes.size(); j++)
                cout << j << ":" << vertex_centralities[nodes[j]] - CB[nodes[j]] << " ";
            cout << endl;
#endif


            for(auto j = CB.begin(); j != CB.end(); j++)
            {
                sensitivityRes[nodes[i]][j->first] = vertex_centralities[j->first] - j->second;
                //nodes[i]->sensitivityValues.push_back(centrality2[j]-vertex_centralities[j]);
            }
#ifdef DEBUG_CHECK_ITERATION
            std::cout<<"finish iteration incremental: " << i << endl;
#endif

            //******
            //std::unordered_map<Node*, std::unordered_map<Node*, std::vector<Node*> > > P_store_;
            for(auto iter = P_store_.begin(); iter != P_store_.end(); iter++)
            {
                for(auto iter_inner = iter->second.begin(); iter_inner != iter->second.end(); iter_inner++)
                {
                    P_store[iter->first][iter_inner->first].clear();
                    for(int ii = 0; ii < iter_inner->second.size(); ii++)
                        P_store[iter->first][iter_inner->first].push_back(iter_inner->second[ii]);
                }
            }
            for(auto iter = P_store_transpose_.begin(); iter != P_store_transpose_.end(); iter++)
            {
                for(auto iter_inner = iter->second.begin(); iter_inner != iter->second.end(); iter_inner++)
                {
                    P_store_transpose[iter->first][iter_inner->first].clear();
                    for(int ii = 0; ii < iter_inner->second.size(); ii++)
                        P_store_transpose[iter->first][iter_inner->first].push_back(iter_inner->second[ii]);
                }
            }
            //std::unordered_map<Node*, std::unordered_map<Node*, int> > sigma_store_;
            for(auto iter = sigma_store_.begin(); iter != sigma_store_.end(); iter++)
            {
                for(auto iter_inner = iter->second.begin(); iter_inner != iter->second.end(); iter_inner++)
                {
                        sigma_store[iter->first][iter_inner->first] = iter_inner->second;
                }
            }
            for(auto iter = sigma_store_transpose_.begin(); iter != sigma_store_transpose_.end(); iter++)
            {
                for(auto iter_inner = iter->second.begin(); iter_inner != iter->second.end(); iter_inner++)
                {
                        sigma_store_transpose[iter->first][iter_inner->first] = iter_inner->second;
                }
            }
            for(auto iter = distance_store_.begin(); iter != distance_store_.end(); iter++)
            {
                for(auto iter_inner = iter->second.begin(); iter_inner != iter->second.end(); iter_inner++)
                {
                        distance_store[iter->first][iter_inner->first] = iter_inner->second;
                }
            }
            for(auto iter = distance_store_transpose_.begin(); iter != distance_store_transpose_.end(); iter++)
            {
                for(auto iter_inner = iter->second.begin(); iter_inner != iter->second.end(); iter_inner++)
                {
                        distance_store_transpose[iter->first][iter_inner->first] = iter_inner->second;
                }
            }
            for(auto iter = cost_store_.begin(); iter != cost_store_.end(); iter++)
            {
                for(auto iter_inner = iter->second.begin(); iter_inner != iter->second.end(); iter_inner++)
                {
                        cost_store[iter->first][iter_inner->first] = iter_inner->second;
                }
            }
            for(auto iter = cost_store_transpose_.begin(); iter != cost_store_transpose_.end(); iter++)
            {
                for(auto iter_inner = iter->second.begin(); iter_inner != iter->second.end(); iter_inner++)
                {
                        cost_store_transpose[iter->first][iter_inner->first] = iter_inner->second;
                }
            }
            for(auto iter = vertex_centralities.begin(); iter != vertex_centralities.end(); iter++)
                CB[iter->first] = iter->second;

            //******


        }
}



void Betweenness::calSensitivityOriginal(std::vector<Node*> &nodes, std::unordered_map<Node*, double> vertex_centralities)
{

    //***********************
    //need to store:
    /*
        std::unordered_map<Node*, std::unordered_map<Node*, std::vector<Node*> > > P_store;
        std::unordered_map<Node*, std::unordered_map<Node*, int> > sigma_store;
        std::unordered_map<Node*, std::unordered_map<Node*, double> > distance_store;
        std::unordered_map<Node*, std::unordered_map<Node*, double> > cost_store;

        std::unordered_map<Node*, std::unordered_map<Node*, std::vector<Node*> > > P_store_transpose;
        std::unordered_map<Node*, std::unordered_map<Node*, int> > sigma_store_transpose;
        std::unordered_map<Node*, std::unordered_map<Node*, double> > distance_store_transpose;
        std::unordered_map<Node*, std::unordered_map<Node*, double> > cost_store_transpose;
     */
    std::unordered_map<Node*, std::unordered_map<Node*, std::vector<Node*> > > P_store_;
    std::unordered_map<Node*, std::unordered_map<Node*, int> > sigma_store_;
    std::unordered_map<Node*, std::unordered_map<Node*, double> > distance_store_;
    std::unordered_map<Node*, std::unordered_map<Node*, double> > cost_store_;

    std::unordered_map<Node*, std::unordered_map<Node*, std::vector<Node*> > > P_store_transpose_;
    std::unordered_map<Node*, std::unordered_map<Node*, int> > sigma_store_transpose_;
    std::unordered_map<Node*, std::unordered_map<Node*, double> > distance_store_transpose_;
    std::unordered_map<Node*, std::unordered_map<Node*, double> > cost_store_transpose_;

    for(auto iter = P_store.begin(); iter != P_store.end(); iter++)
    {
        std::unordered_map<Node*, std::vector<Node*> > P;
        for(auto iter_inner = iter->second.begin(); iter_inner != iter->second.end(); iter_inner++)
        {
            for(int ii = 0; ii < iter_inner->second.size(); ii++)
                P[iter_inner->first].push_back(iter_inner->second[ii]);
        }
        P_store_[iter->first] = P;
    }
    for(auto iter = P_store_transpose.begin(); iter != P_store_transpose.end(); iter++)
    {
        std::unordered_map<Node*, std::vector<Node*> > P;
        for(auto iter_inner = iter->second.begin(); iter_inner != iter->second.end(); iter_inner++)
        {
            for(int ii = 0; ii < iter_inner->second.size(); ii++)
                P[iter_inner->first].push_back(iter_inner->second[ii]);
        }
        P_store_transpose_[iter->first] = P;
    }
    for(auto iter = sigma_store.begin(); iter != sigma_store.end(); iter++)
    {
        std::unordered_map<Node*, int > sigma;
        for(auto iter_inner = iter->second.begin(); iter_inner != iter->second.end(); iter_inner++)
        {
                sigma[iter_inner->first] = iter_inner->second;
        }
        sigma_store_[iter->first] = sigma;
    }
    for(auto iter = sigma_store_transpose.begin(); iter != sigma_store_transpose.end(); iter++)
    {
        std::unordered_map<Node*, int > sigma;
        for(auto iter_inner = iter->second.begin(); iter_inner != iter->second.end(); iter_inner++)
        {
                sigma[iter_inner->first] = iter_inner->second;
        }
        sigma_store_transpose_[iter->first] = sigma;
    }
    for(auto iter = distance_store.begin(); iter != distance_store.end(); iter++)
    {
        std::unordered_map<Node*, double > distance;
        for(auto iter_inner = iter->second.begin(); iter_inner != iter->second.end(); iter_inner++)
        {
                distance[iter_inner->first] = iter_inner->second;
        }
        distance_store_[iter->first] = distance;
    }
    for(auto iter = distance_store_transpose.begin(); iter != distance_store_transpose.end(); iter++)
    {
        std::unordered_map<Node*, double > distance;
        for(auto iter_inner = iter->second.begin(); iter_inner != iter->second.end(); iter_inner++)
        {
                distance[iter_inner->first] = iter_inner->second;
        }
        distance_store_transpose_[iter->first] = distance;
    }
    for(auto iter = cost_store.begin(); iter != cost_store.end(); iter++)
    {
        std::unordered_map<Node*, double > cost;
        for(auto iter_inner = iter->second.begin(); iter_inner != iter->second.end(); iter_inner++)
        {
                cost[iter_inner->first] = iter_inner->second;
        }
        cost_store_[iter->first] = cost;
    }
    for(auto iter = cost_store_transpose.begin(); iter != cost_store_transpose.end(); iter++)
    {
        std::unordered_map<Node*, double > cost;
        for(auto iter_inner = iter->second.begin(); iter_inner != iter->second.end(); iter_inner++)
        {
                cost[iter_inner->first] = iter_inner->second;
        }
        cost_store_transpose_[iter->first] = cost;
    }



    //***********************
        int n = nodes.size();
        for(int i = 0; i < n; i++)
        {
//            brandes_implementation_init(nodes);
//            brandes_implementation_weighted(nodes, false);
//            brandes_implementation_weighted(nodes, true);

            for(int gi = 0; gi < n; gi++)
            {
                float sumWeights = 0;
                for(unsigned int k = 0; k < nodes[gi]->edges.size(); k++)
                {
                    Edge* edge = nodes[gi]->edges[k];
//                    float weight = 1.0;
                    //custom the weight value here!!!
                    sumWeights += edge->weight;
                }
                for(unsigned int k=0;k<nodes[gi]->edges.size();k++) {
                    Edge * edge = nodes[gi]->edges[k];
                    int gj = edge->getNode1()->getIndex();
                    float weight = edge->weight;
                    //custom the weight value here!!!
                    if((i==gi || i==gj) && sumWeights!=0)
                    {
//                        edge->weight = weight - weight/sumWeights;
//                        cout << edge->getNode0()->getIndex() << " " << edge->getNode1()->getIndex() << endl;
//                        cout << weight - weight/sumWeights << endl;
                        //insertEdge(edge->getNode0(), edge->getNode1(), weight - weight/sumWeights);
                        cost_store[edge->getNode0()][edge->getNode1()] = weight - weight/sumWeights;
                    }
                    //edge->weight = ((i==gi || i==gj) && sumWeights!=0)? weight + weight/sumWeights: weight;
                }
            }

            brandes_implementation_weighted(nodes, false);

#ifdef DEBUG_SENSITIVITY
            cout << "sensitivity original: " << i << endl;
//            for(auto j = CB.begin(); j != CB.end(); j++)
//            {
//                cout << j->first->getIndex() << ":" << vertex_centralities[j->first] - j->second << " ";
//            }
            for(int j = 0; j < nodes.size(); j++)
                cout << j << ":" << vertex_centralities[nodes[j]] - CB[nodes[j]] << " ";
            cout << endl;
#endif


            for(auto j = CB.begin(); j != CB.end(); j++)
            {
                sensitivityRes[nodes[i]][j->first] = vertex_centralities[j->first] - j->second;
                //nodes[i]->sensitivityValues.push_back(centrality2[j]-vertex_centralities[j]);
            }
#ifdef DEBUG_CHECK_ITERATION
            std::cout<<"finish iteration incremental: " << i << endl;
#endif
            //******
            //std::unordered_map<Node*, std::unordered_map<Node*, std::vector<Node*> > > P_store_;
            for(auto iter = P_store_.begin(); iter != P_store_.end(); iter++)
            {
                for(auto iter_inner = iter->second.begin(); iter_inner != iter->second.end(); iter_inner++)
                {
                    P_store[iter->first][iter_inner->first].clear();
                    for(int ii = 0; ii < iter_inner->second.size(); ii++)
                        P_store[iter->first][iter_inner->first].push_back(iter_inner->second[ii]);
                }
            }
            for(auto iter = P_store_transpose_.begin(); iter != P_store_transpose_.end(); iter++)
            {
                for(auto iter_inner = iter->second.begin(); iter_inner != iter->second.end(); iter_inner++)
                {
                    P_store_transpose[iter->first][iter_inner->first].clear();
                    for(int ii = 0; ii < iter_inner->second.size(); ii++)
                        P_store_transpose[iter->first][iter_inner->first].push_back(iter_inner->second[ii]);
                }
            }
            //std::unordered_map<Node*, std::unordered_map<Node*, int> > sigma_store_;
            for(auto iter = sigma_store_.begin(); iter != sigma_store_.end(); iter++)
            {
                for(auto iter_inner = iter->second.begin(); iter_inner != iter->second.end(); iter_inner++)
                {
                        sigma_store[iter->first][iter_inner->first] = iter_inner->second;
                }
            }
            for(auto iter = sigma_store_transpose_.begin(); iter != sigma_store_transpose_.end(); iter++)
            {
                for(auto iter_inner = iter->second.begin(); iter_inner != iter->second.end(); iter_inner++)
                {
                        sigma_store_transpose[iter->first][iter_inner->first] = iter_inner->second;
                }
            }
            for(auto iter = distance_store_.begin(); iter != distance_store_.end(); iter++)
            {
                for(auto iter_inner = iter->second.begin(); iter_inner != iter->second.end(); iter_inner++)
                {
                        distance_store[iter->first][iter_inner->first] = iter_inner->second;
                }
            }
            for(auto iter = distance_store_transpose_.begin(); iter != distance_store_transpose_.end(); iter++)
            {
                for(auto iter_inner = iter->second.begin(); iter_inner != iter->second.end(); iter_inner++)
                {
                        distance_store_transpose[iter->first][iter_inner->first] = iter_inner->second;
                }
            }
            for(auto iter = cost_store_.begin(); iter != cost_store_.end(); iter++)
            {
                for(auto iter_inner = iter->second.begin(); iter_inner != iter->second.end(); iter_inner++)
                {
                        cost_store[iter->first][iter_inner->first] = iter_inner->second;
                }
            }
            for(auto iter = cost_store_transpose_.begin(); iter != cost_store_transpose_.end(); iter++)
            {
                for(auto iter_inner = iter->second.begin(); iter_inner != iter->second.end(); iter_inner++)
                {
                        cost_store_transpose[iter->first][iter_inner->first] = iter_inner->second;
                }
            }
            for(auto iter = vertex_centralities.begin(); iter != vertex_centralities.end(); iter++)
                CB[iter->first] = iter->second;

            //******
        }
}


void Betweenness::insertEdge(Node* src, Node* dest, double cost)
{

//    std::map<Node*, double> CB;
//
//    std::map<Node*, std::map<Node*, std::vector<Node*> > > P_store;
//    std::map<Node*, std::map<Node*, int> > sigma_store;
//    std::map<Node*, std::map<Node*, double> > distance_store;
//    std::map<Node*, std::map<Node*, double> > cost_store;

#ifdef DEBGU_GET_STORE_VALUES
    cout << "CB::" << endl;
    for(auto i = CB.begin(); i != CB.end(); i++)
    {
        cout << i->first->getIndex() << ":" << i->second << endl;
    }
    cout << "CB_transpose::" << endl;
    for(auto i = CB_transpose.begin(); i != CB_transpose.end(); i++)
    {
        cout << i->first->getIndex() << ":" << i->second << endl;
    }
    cout << "P::" << endl;
    for(auto i = P_store.begin(); i != P_store.end(); i++)
    {
        for(auto j = i->second.begin(); j != i->second.end(); j++)
        {
            cout << "<" << i->first->getIndex() << "," << j->first->getIndex() << ">: ";
            for(int k = 0; k < j->second.size(); k++)
                cout << j->second[k]->getIndex() << " ";
            cout << endl;
        }
    }
    cout << "P_transpose::" << endl;
    for(auto i = P_store_transpose.begin(); i != P_store_transpose.end(); i++)
    {
        for(auto j = i->second.begin(); j != i->second.end(); j++)
        {
            cout << "<" << i->first->getIndex() << "," << j->first->getIndex() << ">: ";
            for(int k = 0; k < j->second.size(); k++)
                cout << j->second[k]->getIndex() << " ";
            cout << endl;
        }
    }

    cout << "sigma::" << endl;
    for(auto i = sigma_store.begin(); i != sigma_store.end(); i++)
    {
        for(auto j = i->second.begin(); j != i->second.end(); j++)
        {
            cout << "<" << i->first->getIndex() << "," << j->first->getIndex() << ">: " << j->second << endl;
        }
    }
    cout << "sigma_transpose::" << endl;
    for(auto i = sigma_store_transpose.begin(); i != sigma_store_transpose.end(); i++)
    {
        for(auto j = i->second.begin(); j != i->second.end(); j++)
        {
            cout << "<" << i->first->getIndex() << "," << j->first->getIndex() << ">: " << j->second << endl;
        }
    }

    cout << "distance::" << endl;
    for(auto i = distance_store.begin(); i != distance_store.end(); i++)
    {
        for(auto j = i->second.begin(); j != i->second.end(); j++)
        {
            cout << "<" << i->first->getIndex() << "," << j->first->getIndex() << ">: " << j->second << endl;
        }
    }
    cout << "distance_transpose::" << endl;
    for(auto i = distance_store_transpose.begin(); i != distance_store_transpose.end(); i++)
    {
        for(auto j = i->second.begin(); j != i->second.end(); j++)
        {
            cout << "<" << i->first->getIndex() << "," << j->first->getIndex() << ">: " << j->second << endl;
        }
    }

    cout << "cost::" << endl;
    for(auto i = cost_store.begin(); i != cost_store.end(); i++)
    {
        for(auto j = i->second.begin(); j != i->second.end(); j++)
        {
            cout << "<" << i->first->getIndex() << "," << j->first->getIndex() << ">: " << j->second << endl;
        }
    }
    cout << "cost_transpose::" << endl;
    for(auto i = cost_store_transpose.begin(); i != cost_store_transpose.end(); i++)
    {
        for(auto j = i->second.begin(); j != i->second.end(); j++)
        {
            cout << "<" << i->first->getIndex() << "," << j->first->getIndex() << ">: " << j->second << endl;
        }
    }
#endif


    sigma_old.clear();
    distance_old.clear();
    trackLost.clear();
    pairsDone.clear();
    sigma_old_transpose.clear();
    distance_old_transpose.clear();
    trackLost_transpose.clear();
    pairsDone_transpose.clear();

    //mind this, undirected graph...
    cost_store[src][dest] = cost;
    cost_store_transpose[dest][src] = cost;

//    cost_store[dest][src] = cost;

//    TimeLogger* logger = TimeLogger::Instance();
//    logger->start();

    std::vector<Node*> sinks = insertUpdate(dest, src, src, true);
    std::vector<Node*> sources = insertUpdate(src, dest, dest, false);
//    logger->markIt("after computing source and target: ");

    for(int i = 0; i < sinks.size(); i++)
        insertUpdate(src, dest, sinks[i], false);
//    logger->markIt("after finish target: ");
    for(int i = 0; i < sources.size(); i++)
        insertUpdate(dest, src, sources[i], true);
//    logger->markIt("after finish source: ");

//    increaseBetweenness(true);
    increaseBetweenness(false);
//    logger->markIt("after finish increase betweenness: ");
//    logger->outputToScreen();
}
std::vector<Node*> Betweenness::insertUpdate(Node* src, Node* dest, Node* z, bool isTranspose)
{
    std::vector<std::pair<Node*, Node*> > workSet;
    std::vector<Node*> visitedVertices;
    std::vector<Node*> affectedVertices;

//    std::cout << "cost from src to dest::" << cost_store[src][dest] << std::endl;
//    std::cout << "distance from src to dest::" << distance_store[src][dest] << std::endl;

    std::unordered_map<Node*, std::unordered_map<Node*, std::vector<Node*> > > &P_now = isTranspose? P_store_transpose:P_store;
    std::unordered_map<Node*, std::unordered_map<Node*, int> > &sigma_now = isTranspose? sigma_store_transpose:sigma_store;
    std::unordered_map<Node*, std::unordered_map<Node*, double> > &distance_now = isTranspose? distance_store_transpose:distance_store;
    std::unordered_map<Node*, std::unordered_map<Node*, double> > &cost_now = isTranspose? cost_store_transpose:cost_store;

    std::unordered_map<Node*, std::unordered_map<Node*, int> > &sigma_old_now = isTranspose? sigma_old_transpose:sigma_old;
    std::unordered_map<Node*, std::unordered_map<Node*, double> > &distance_old_now = isTranspose? distance_old_transpose:distance_old;
    std::vector<std::pair<Node*, Node*> > &pairsDone_now = /*isTranspose? pairsDone_transpose :*/ pairsDone;
    std::vector<std::tuple<Node*, Node*, Node*> > &trackLost_now = isTranspose? trackLost_transpose:trackLost;

    workSet.push_back(make_pair(src, dest));
    visitedVertices.push_back(src);

    while (workSet.size() != 0) {   // 4
        Node* x = workSet.back().first;
        Node* y = workSet.back().second;
        workSet.pop_back();
        double alt = cost_now[x][y] + distance_now[y][z];

        if(alt < distance_now[x][z])
        {
            if(!isIn(make_pair(x, z), sigma_old_now))
            {
                distance_old_now[x][z] = distance_now[x][z];
                sigma_old_now[x][z] = sigma_now[x][z];
                reduceBetweenness(x, z, false);    // 10
                sigma_now[x][z] = 0;
                P_now[x][z].clear();
            }
            if(isIn(make_pair(x, z), pairsDone_now))
                pairsDone_now.erase(find(pairsDone_now.begin(), pairsDone_now.end(), make_pair(x, z)));
            distance_now[x][z] = alt;
        }
        if(alt == distance_now[x][z] && distance_now[x][z] != DBL_MAX)
        {
            if(!isIn(make_pair(x, z), pairsDone_now))
            {
                if(!isIn(make_pair(x, z), sigma_old_now))
                    reduceBetweenness(x, z, false);    // 18
                if(sigma_now[x][z] != 0)
                    sigma_old_now[x][z] = sigma_now[x][z];

                sigma_now[x][z] = sigma_now[x][z] + (sigma_now[x][src] * 1 * sigma_now[dest][z]);
                if(std::find(P_now[x][y].begin(), P_now[x][y].end(), x) == P_now[x][y].end())
                    P_now[x][y].push_back(x);
                for(int iter = 0; iter < P_now[y][z].size(); iter++)
                {
                    if(std::find(P_now[x][z].begin(), P_now[x][z].end(), P_now[y][z][iter]) == P_now[x][z].end())
                        P_now[x][z].push_back(P_now[y][z][iter]);
                }
            }
            pairsDone_now.push_back(make_pair(x, z));
            affectedVertices.push_back(x);
            // need pred(x)...
            auto iter = isTranspose?x->pred_transpose.begin():x->pred.begin();
            auto iter_end = isTranspose?x->pred_transpose.end():x->pred.end();
            for(; iter != iter_end; iter++)
            {
                auto pred_vec = iter->second;
                for (int inner = 0; inner < pred_vec.size(); inner++) {
                    Node* u = pred_vec[inner];
                    if(SP(u, x, src, isTranspose) && std::find(visitedVertices.begin(), visitedVertices.end(), u) == visitedVertices.end())
                    {
                        workSet.push_back(make_pair(u, x));
                        visitedVertices.push_back(u);
                    }
                }

            }
        }
    }
    return affectedVertices;
}
double Betweenness::getDistanceOldVal(Node *x, Node *y, bool isTranspose)
{
    if(isIn(std::make_pair(x, y), isTranspose?distance_old_transpose:distance_old))
        return isTranspose?distance_old_transpose[x][y]:distance_old[x][y];
    else if(isTranspose)
        return distance_store_transpose[x][y];
    return distance_store[x][y];
}
int Betweenness::getSigmaOldVal(Node *x, Node *y, bool isTranspose)
{
    if(isIn(std::make_pair(x, y), isTranspose? sigma_old_transpose:sigma_old))
        return isTranspose? sigma_old_transpose[x][y]:sigma_old[x][y];
    else if(isTranspose)
        return sigma_store_transpose[x][y];
    return sigma_store[x][y];
}

void Betweenness::reduceBetweenness(Node* a, Node* z, bool isTranspose)
{

    std::unordered_map<Node*, std::unordered_map<Node*, std::vector<Node*> > > &P_now = isTranspose? P_store_transpose:P_store;
    std::unordered_map<Node*, std::unordered_map<Node*, int> > &sigma_now = isTranspose? sigma_store_transpose:sigma_store;
    std::unordered_map<Node*, std::unordered_map<Node*, double> > &distance_now = isTranspose? distance_store_transpose:distance_store;
    std::unordered_map<Node*, std::unordered_map<Node*, double> > &cost_now = isTranspose? cost_store_transpose:cost_store;

    std::unordered_map<Node*, std::unordered_map<Node*, int> > &sigma_old_now = isTranspose? sigma_old_transpose:sigma_old;
    std::unordered_map<Node*, std::unordered_map<Node*, double> > &distance_old_now = isTranspose? distance_old_transpose:distance_old;
    std::vector<std::pair<Node*, Node*> > &pairsDone_now = /*isTranspose? pairsDone_transpose :*/ pairsDone;
    std::vector<std::tuple<Node*, Node*, Node*> > &trackLost_now = isTranspose? trackLost_transpose:trackLost;


    if(getSigmaOldVal(a, z, isTranspose) == 0)
        return;
    std::vector<Node*> known;
    std::vector<Node*> stack;

//    if(a->getIndex() == 0 && z->getIndex() == 3)
//        std::cout << "here" << endl;

    for(int i = 0; i < P_now[a][z].size(); i++)
    {
        Node* n = P_now[a][z][i];
        if(distance_now[a][z] != getDistanceOldVal(a, n, isTranspose) + getDistanceOldVal(n, z, isTranspose))
            continue;
        else if(a != n && n != z)
        {
//            cout << "node" << n->getIndex() << endl;
//            cout << "line num:: 8" << endl;
//            cout << "from " << a->getIndex() << " to " << z->getIndex() << endl;
//            std::cout << getSigmaOldVal(a, n, isTranspose) * getSigmaOldVal(n, z, isTranspose) << endl;
//            std::cout << getSigmaOldVal(a, z, isTranspose) << endl;
//            std::cout << CB[n] << endl;
            CB[n] = CB[n] - ((getSigmaOldVal(a, n, isTranspose) * getSigmaOldVal(n, z, isTranspose)*1.0) / getSigmaOldVal(a, z, isTranspose));  //8
//            std::cout << CB[n] << endl;
//            cout << "8:: values of <a,z,n>:" << a->getIndex() << ", " << z->getIndex() << ", " << n->getIndex() << endl;
            trackLost_now.push_back(make_tuple(a, z, n));
        }
        stack.push_back(n);
        known.push_back(n);
    }
    while (stack.size() != 0) {
        Node* p = stack.back();
        stack.pop_back();
        known.push_back(p);
        for(int i = 0; i < P_now[a][p].size(); i++)
        {
            Node* n = P_now[a][p][i];
            if(distance_now[a][z] != getDistanceOldVal(a, n, isTranspose) + getDistanceOldVal(n, z, isTranspose))
                continue;
            else if(a != n && n != z && std::find(known.begin(), known.end(), n) == known.end())
            {
//                cout << "node" << n->getIndex() << endl;
//                cout << "line num:: 18" << endl;
//                cout << "from " << a->getIndex() << " to " << z->getIndex() << endl;
//                std::cout << getSigmaOldVal(a, n, isTranspose) * getSigmaOldVal(n, z, isTranspose) << endl;
//                std::cout << getSigmaOldVal(a, z, isTranspose) << endl;
//                std::cout << CB[n] << endl;
                CB[n] = CB[n] - ((getSigmaOldVal(a, n, isTranspose) * getSigmaOldVal(n, z, isTranspose) * 1.0) / getSigmaOldVal(a, z, isTranspose));  //18
//                std::cout << CB[n] << endl;
//                cout << "18:: values of <a,z,n>:" << a->getIndex() << ", " << z->getIndex() << ", " << n->getIndex() << endl;
                trackLost_now.push_back(make_tuple(a, z, n));
                stack.push_back(n);
                known.push_back(n);
            }
        }
    }
    std::vector<Node*> alreadyKnown;
    for(int i = 0; i < known.size(); i++)
        alreadyKnown.push_back(known[i]);
    alreadyKnown.push_back(a);
//    alreadyKnown.push_back(z);



    for(int i = 0; i < trackLost_now.size(); i++)
    {
        Node* v = std::get<2>(trackLost_now[i]);
        Node* v1 = std::get<0>(trackLost_now[i]);
        Node* v2 = std::get<1>(trackLost_now[i]);

//        if(std::find(known.begin(), known.end(), v1) != known.end() &&
//           std::find(known.begin(), known.end(), v2) != known.end() &&
//           distance_now[a][z] == getDistanceOldVal(a, v, isTranspose) + getDistanceOldVal(v, z, isTranspose))
        if(std::find(alreadyKnown.begin(), alreadyKnown.end(), v1) != alreadyKnown.end() &&
           std::find(alreadyKnown.begin(), alreadyKnown.end(), v2) != alreadyKnown.end() &&
           distance_now[a][z] == getDistanceOldVal(a, v, isTranspose) + getDistanceOldVal(v, z, isTranspose))
        {
            if(std::find(alreadyKnown.begin(), alreadyKnown.end(), v) == alreadyKnown.end())
            {
//                cout << "node" << v->getIndex() << endl;
//                cout << "line num:: 24" << endl;
//                cout << "from " << a->getIndex() << " to " << z->getIndex() << endl;
//                std::cout << getSigmaOldVal(a, v, isTranspose) * getSigmaOldVal(v, z, isTranspose) << endl;
//                std::cout << getSigmaOldVal(a, z, isTranspose) << endl;
//                std::cout << CB[v] << endl;
                CB[v] = CB[v] - ((getSigmaOldVal(a, v, isTranspose) * getSigmaOldVal(v, z, isTranspose) * 1.0) / getSigmaOldVal(a, z, isTranspose));  //24
//                std::cout << CB[v] << endl;

                alreadyKnown.push_back(v);
//                cout << "24::values of <a,z,n>:" << a->getIndex() << ", " << z->getIndex() << ", " << v->getIndex() << endl;
//                cout << "v1 and v2::" << v1->getIndex() << ", " << v2->getIndex() << endl;
                trackLost_now.push_back(std::make_tuple(a, z, v));
            }
        }
    }

}
void Betweenness::increaseBetweenness(bool isTranspose)
{
//    for(int i = 0; i < trackLost.size(); i++)
//    {
//        Node* n = std::get<2>(trackLost[i]);
//        Node* a = std::get<0>(trackLost[i]);
//        Node* z = std::get<1>(trackLost[i]);
//        cout << "<" << a->getIndex() << "," << z->getIndex() << "," << n->getIndex() << ">" << endl;
//    }

//    std::cout << "****** begin increase ******" << endl;

    std::unordered_map<Node*, std::unordered_map<Node*, std::vector<Node*> > > &P_now = isTranspose? P_store_transpose:P_store;
    std::unordered_map<Node*, std::unordered_map<Node*, int> > &sigma_now = isTranspose? sigma_store_transpose:sigma_store;
    std::unordered_map<Node*, std::unordered_map<Node*, double> > &distance_now = isTranspose? distance_store_transpose:distance_store;
    std::unordered_map<Node*, std::unordered_map<Node*, double> > &cost_now = isTranspose? cost_store_transpose:cost_store;

    std::unordered_map<Node*, std::unordered_map<Node*, int> > &sigma_old_now = isTranspose? sigma_old_transpose:sigma_old;
    std::unordered_map<Node*, std::unordered_map<Node*, double> > &distance_old_now = isTranspose? distance_old_transpose:distance_old;
    std::vector<std::pair<Node*, Node*> > &pairsDone_now = /*isTranspose? pairsDone_transpose :*/ pairsDone;
    std::vector<std::tuple<Node*, Node*, Node*> > &trackLost_now = isTranspose? trackLost_transpose:trackLost;

    for(auto outer = sigma_old_now.begin(); outer != sigma_old_now.end(); outer++)
    {
        auto src = outer->first;
        for(auto inner = outer->second.begin(); inner != outer->second.end(); inner++)
        {
            auto dest = inner->first;
            std::vector<Node*> known;
            std::vector<Node*> stack;
//            std::cout << "******" << endl;
//            cout << "P::" << src->getIndex() <<" to "<< dest->getIndex() << endl;
//            for(int i = 0; i < P_store[src][dest].size(); i++)
//                cout << P_store[src][dest][i]->getIndex() << " ";
//            cout << endl;
            for(int i = 0; i < P_now[src][dest].size(); i++)
            {
                Node* n = P_now[src][dest][i];

                stack.push_back(n);
                known.push_back(n);
                if(src != n && n != dest)
                {
//                    cout << "node" << n->getIndex() << endl;
//                    cout << "path::" << src->getIndex() << "->" << dest->getIndex() << endl;
//                    cout << "line num:: 6" << endl;
//                    std::cout << sigma_store[src][n] * sigma_store[n][dest] << endl;
//                    std::cout << sigma_store[src][dest] << endl;
//                    std::cout << CB[n] << endl;

                    CB[n] = CB[n] + (((sigma_now[src][n] * sigma_now[n][dest]) * 1.0) / sigma_now[src][dest]); // 6
//                    std::cout << CB[n];
//                    std::cout << endl;

                }

            }
            while (stack.size() != 0) {
                Node* n = stack.back();
                stack.pop_back();
                known.push_back(n);


//                std::cout << "******" << endl;
//                cout << "P::" << src->getIndex() <<" to "<< dest->getIndex() << endl;
//                for(int i = 0; i < P_store[src][dest].size(); i++)
//                    cout << P_store[src][dest][i]->getIndex() << " ";
//                cout << endl;

                for(int i = 0; i < P_now[src][n].size(); i++)
                {
                    Node* p = P_now[src][n][i];
                    if(p != src && p != dest && std::find(known.begin(), known.end(), p) == known.end())
                    {
                        stack.push_back(p);
                        known.push_back(p);

//                        cout << "node" << p->getIndex() << endl;
//                        cout << "path::" << src->getIndex() << "->" << dest->getIndex() << endl;
//                        cout << "line num:: 13" << endl;
//                        std::cout << sigma_store[src][p] * sigma_store[p][dest] << endl;
//                        std::cout << sigma_store[src][dest] << endl;
//                        std::cout << CB[p] << endl;

                        CB[p] = CB[p] + ((sigma_now[src][p] * sigma_now[p][dest])*1.0 / sigma_now[src][dest]);  //13

//                        std::cout << CB[p];
//                        std::cout << endl;
                    }
                }
            }
        }
    }
}


//void Betweenness::deleteEdge(Node *src, Node *dest)
//{
//    sigma_old.clear();
//    distance_old.clear();
//    trackLost.clear();
//    pairsDone.clear();

////    cost_store[src][dest] = DBL_MAX;

//    std::vector<Node*> sinks = deleteUpdate(dest, src, src);
//    std::vector<Node*> sources = deleteUpdate(src, dest, dest);
//    for(int i = 0; i < sinks.size(); i++)
//        deleteUpdate(src, dest, sinks[i]);
//    for(int i = 0; i < sources.size(); i++)
//        deleteUpdate(dest, src, sources[i]);
//    increaseBetweenness();
//}

//std::vector<Node*> Betweenness::deleteUpdate(Node* src, Node* dest, Node* z)
//{
//    std::vector<std::pair<Node*, Node*> > workSet;
//    std::vector<Node*> visitedVertices;
//    std::vector<Node*> affectedVertices;

////    std::cout << "cost from src to dest::" << cost_store[src][dest] << std::endl;
////    std::cout << "distance from src to dest::" << distance_store[src][dest] << std::endl;

//    workSet.push_back(make_pair(src, dest));
//    visitedVertices.push_back(src);

//    while (workSet.size() != 0) {
//        Node* x = workSet.back().first;
//        Node* y = workSet.back().second;
//        workSet.pop_back();
//        double alt = cost_store[x][y] + distance_store[y][z];
//        if(alt > distance_store[x][z])
//        {
//            if(!isIn(make_pair(x, z), sigma_old))
//            {
//                distance_old[x][z] = distance_store[x][z];
//                sigma_old[x][z] = sigma_store[x][z];
//                reduceBetweenness(x, z);
//                sigma_store[x][z] = 0;
//                P_store[x][z].clear();

//                std::map<double, vector<Node*> > tmp_store;
//                for(int iter = 0; iter < x->edges.size(); iter++)
//                {
//                    Node* v = x->edges[iter]->getNode1();
//                    if(cost_store[x][v] != DBL_MAX && distance_store[v][z] != DBL_MAX)
//                        tmp_store[distance_store[v][z] + cost_store[x][v]].push_back(v);
//                }
//                if(tmp_store.size() != 0)
//                {
//                    for(int iter = 0; iter < tmp_store[0].size(); iter++)
//                    {
//                        Node* v = tmp_store[0][iter];
//                        sigma_store[x][z] += sigma_store[v][z];
//                        for(int p_iter = 0; p_iter < P_store[v][z].size(); p_iter++)
//                            P_store[x][z].push_back(P_store[v][z][p_iter]);
//                    }
//                    alt = tmp_store.begin()->first;
//                }

//            }
//            if(isIn(make_pair(x, z), pairsDone))
//                pairsDone.erase(find(pairsDone.begin(), pairsDone.end(), make_pair(x, z)));
//            distance_store[x][z] = alt;
//        }
//        if(alt == distance_store[x][z] && distance_store[x][z] != DBL_MAX)
//        {
//            cost_store[src][dest] = DBL_MAX;
//            if(!isIn(make_pair(x, z), pairsDone))
//            {
//                if(!isIn(make_pair(x, z), sigma_old))
//                {
//                    reduceBetweenness(x, z);
//                    //if(sigma_store[x][z] != 0)
//                    sigma_old[x][z] = sigma_store[x][z];
//                    sigma_store[x][z] = sigma_store[x][z] - (sigma_store[x][src] * 1 * sigma_store[dest][z]);
//                    if(sigma_store[x][z] == 0)
//                    {
//                        P_store[x][z].clear();

//                        std::map<double, vector<Node*> > tmp_store;
//                        for(int iter = 0; iter < x->edges.size(); iter++)
//                        {
//                            Node* v = x->edges[iter]->getNode1();
//                            if(cost_store[x][v] != DBL_MAX && distance_store[v][z] != DBL_MAX)
//                                tmp_store[distance_store[v][z] + cost_store[x][v]].push_back(v);
//                        }
//                        if(tmp_store.size() != 0)
//                        {
//                            for(int iter = 0; iter < tmp_store[0].size(); iter++)
//                            {
//                                Node* v = tmp_store[0][iter];
//                                sigma_store[x][z] += sigma_store[v][z];
//                                for(int p_iter = 0; p_iter < P_store[v][z].size(); p_iter++)
//                                    P_store[x][z].push_back(P_store[v][z][p_iter]);
//                            }
//                            distance_store[x][z] = tmp_store.begin()->first;
//                        }
//                        else
//                            distance_store[x][z] = DBL_MAX;
//                    }
//                    else //remove P[y][z] from P[x][z]
//                    {
//                        for(int iter = 0; iter < P_store[y][z].size(); iter++)
//                        {
//                            P_store[x][z].erase(find(P_store[y][z].begin(), P_store[y][z].end(), P_store[y][z][iter]));
//                        }
//                    }
//                }

//            }
//            pairsDone.push_back(make_pair(x, z));
//            affectedVertices.push_back(x);
//            // need pred(x)...
//            for(auto iter = x->pred.begin(); iter != x->pred.end(); iter++)
//            {
//                auto pred_vec = iter->second;
//                for (int inner = 0; inner < pred_vec.size(); inner++) {
//                    Node* u = pred_vec[inner];
//                    if(SP(u, x, src) && std::find(visitedVertices.begin(), visitedVertices.end(), u) == visitedVertices.end())
//                    {
//                        workSet.push_back(make_pair(u, x));
//                        visitedVertices.push_back(u);
//                    }
//                }

//            }
//        }
//    }
//    return affectedVertices;
//}

bool Betweenness::isIn(std::pair<Node*, Node*> key, std::unordered_map<Node*, std::unordered_map<Node*, double> > &container)
{
    if(container.find(key.first) != container.end())
        if(container[key.first].find(key.second) != container[key.first].end())
            return true;
    return false;
}
bool Betweenness::isIn(std::pair<Node*, Node*> key, std::unordered_map<Node*, std::unordered_map<Node*, int> > &container)
{
    if(container.find(key.first) != container.end())
        if(container[key.first].find(key.second) != container[key.first].end())
            return true;
    return false;
}
bool Betweenness::isIn(std::pair<Node*, Node*> value, std::vector<std::pair<Node*, Node*> > &container)
{
    if(find(container.begin(), container.end(), value) != container.end())
        return true;
    return false;
}
bool Betweenness::isIn(std::tuple<Node*, Node*, Node*> value, std::vector<std::tuple<Node*, Node*, Node*> > &container)
{
    if(find(container.begin(), container.end(), value) != container.end())
        return true;
    return false;
}
bool Betweenness::SP(Node *x, Node *y, Node *z, bool isTranspose)
{
    if(!isTranspose)
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
    else
    {
        if(isIn(make_pair(x, z), distance_store_transpose))
        {
            if(distance_store_transpose[x][z] != DBL_MAX)
            {
                if(distance_store_transpose[x][z] == cost_store_transpose[x][y] + distance_store_transpose[y][z])
                    return true;
            }
        }
        return false;
    }

}

