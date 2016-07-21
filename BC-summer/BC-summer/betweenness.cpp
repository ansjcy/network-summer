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
                                   .vertex_index_map(get(vertex_index, g)).weight_map(get(edge_weight, g)));
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









































