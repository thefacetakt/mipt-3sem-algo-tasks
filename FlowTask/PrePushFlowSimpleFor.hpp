#ifndef _PRE_PUSH_FLOW_SIMPLE_FOR_
#define _PRE_PUSH_FLOW_SIMPLE_FOR_

#include "Edges.hpp"
#include "Net.hpp"
#include "MaxFlowFinder.hpp"
#include <vector>

class PrePushFlowSimpleFor : public MaxFlowFinder {
    Net net_;
    
    std::vector<std::vector<InnerNetEdge>::iterator> firstUnsaturatedEdge_;
    std::vector<unsigned int> height_;
    std::vector<unsigned long long> excess_;

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