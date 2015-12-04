#ifndef _MAX_FLOW_FINDER_
#define _MAX_FLOW_FINDER_

#include "MaxFlowDescription.hpp"
#include "Net.hpp"
#include <vector>

class MaxFlowFinder {
protected:
    Net net_;
    virtual void init(unsigned int verticesCount) = 0;
    virtual void findMaxFlow() = 0;
    virtual void cleanUp() = 0;
public:
    MaxFlowDescription run(unsigned int verticesCount, const std::vector<DirectedEdgeWithStart> &net, unsigned int source, unsigned int sink);
    virtual ~MaxFlowFinder() {
    }
};



#endif