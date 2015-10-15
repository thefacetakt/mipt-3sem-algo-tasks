#ifndef _MAX_FLOW_FINDER_
#define _MAX_FLOW_FINDER_

#include "MaxFlowDescription.hpp"
#include <vector>

class MaxFlowFinder {
public:
    virtual MaxFlowDescription findMaxFlow(unsigned int verticeCount, const std::vector<DirectedEdgeWithStart> &net, unsigned int source, unsigned int sink) = 0;
    virtual ~MaxFlowFinder() {
    }
};



#endif