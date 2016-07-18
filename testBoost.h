///**************************************************************************
// NetZen :  Network Sensitivity Analysis and Visualization
// -------------------
// copyright            : (C) 2009 by Carlos D. Correa
// email                : correac@cs.ucdavis.edu
// ***************************************************************************/
//
///***************************************************************************
// *                                                                         *
// * Copyright: All rights reserved.  May not be used, modified, or copied   *
// * without permission.                                                     *
// *                                                                         *
// ***************************************************************************/
//#ifndef CENTRALITY_H
//#define CENTRALITY_H
//
//#include "object.h"
//#include "analysisTool.h"
//#include "graph.h"
//#include "controller.h"
//#include <vector>
//using namespace std;
//class Node;
//
//#ifdef USINGLAPACK
//class EigenvectorCentrality: public AnalysisTool {
//public:
//    EigenvectorCentrality(Controller *controller);
//    int compute(vector<Node*> &nodes, bool needDerivs);
//    virtual void compute();
//};
//
//class MarkovCentrality: public AnalysisTool {
//public:
//    MarkovCentrality(Controller *controller);
//    int compute(vector<Node*> &nodes, bool needDerivs);
//    virtual void compute();
//};
//#endif //USINGLAPACK
//
//class Betweenness: public AnalysisTool {
//public:
//    Betweenness(Controller *controller);
//    int compute(vector<Node*> &nodes, bool needDerivs);
//    virtual void compute();
//};
//
//class Closeness: public AnalysisTool {
//public:
//    Closeness(Controller *controller);
//    int compute(vector<Node*> &nodes, bool needDerivs);
//    virtual void compute();
//};
//
//class EigenvectorCentralityEigen: public AnalysisTool {
//public:
//    EigenvectorCentralityEigen(Controller *controller);
//    int compute(vector<Node*> &nodes, bool needDerivs);
//    virtual void compute();
//};
//
//class MarkovCentralityEigen: public AnalysisTool {
//public:
//    MarkovCentralityEigen(Controller *controller);
//    int compute(vector<Node*> &nodes, bool needDerivs);
//    virtual void compute();
//};
//
//#endif
