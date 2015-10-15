#ifndef _PRE_PUSH_FLOW_SIMPLE_FOR_
#define _PRE_PUSH_FLOW_SIMPLE_FOR_

#include "Edges.hpp"
#include "MaxFlowFinder.hpp"
#include <vector>

class PrePushFlowSimpleFor : public MaxFlowFinder {
    std::vector<std::vector<InnerNetEdge> > graph;
    std::vector<std::vector<InnerNetEdge>::iterator> firstUnsaturatedEdge;
    std::vector<unsigned int> height;
    std::vector<unsigned long long> excess;
    unsigned int source;
    unsigned int sink;
    unsigned int verticeCount;
    unsigned int edgesCount;
    
    void discharge(unsigned int v);
    
    void relabel(unsigned int v);
    
    void push(unsigned int v, InnerNetEdge &e);
    
    void cleanUp();
    
    MaxFlowDescription findMaxFlowInitialised();
    
public:
    PrePushFlowSimpleFor();
    
    MaxFlowDescription findMaxFlow(unsigned int verticeCount, const std::vector<DirectedEdgeWithStart> &net, unsigned int source, unsigned int sink);
    
    virtual ~PrePushFlowSimpleFor() {
    }
};

#endif