//#include <iostream>
//#include "betweenness.hpp"
//using namespace std;
//int main()
//{
//    Betweenness bc;
//    vector<Node*> nodes;
//    vector<Edge*> edges;
//    for(int i = 0; i < 9; i++)
//    {
//        Node* tmpNode = new Node();
//        tmpNode->setIndex(i);
//        nodes.push_back(tmpNode);
//    }
//    
//    // should add all the edges on the graph!!!!
//    nodes[0]->addEdge(nodes[1]);
//    nodes[0]->addEdge(nodes[2]);
//    nodes[0]->addEdge(nodes[3]);
//    nodes[0]->addEdge(nodes[4]);
//    
//    nodes[1]->addEdge(nodes[0]);
//    
//    nodes[2]->addEdge(nodes[0]);
//    
//    nodes[3]->addEdge(nodes[0]);
//    
//    nodes[4]->addEdge(nodes[0]);
//    nodes[4]->addEdge(nodes[5]);
//    nodes[4]->addEdge(nodes[6]);
//    nodes[4]->addEdge(nodes[7]);
//    
//    nodes[5]->addEdge(nodes[4]);
//    nodes[6]->addEdge(nodes[4]);
//    
//    nodes[7]->addEdge(nodes[4]);
//    nodes[7]->addEdge(nodes[8]);
//    
//    nodes[8]->addEdge(nodes[7]);
//    
//    bc.compute(nodes, true);
//    
////    for (int i = 0; i < 7; i++) {
////        cout<<nodes[i]->index<<endl;
////    }
//    
//}