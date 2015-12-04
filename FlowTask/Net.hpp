#ifndef _NET
#define _NET

#include <vector>
#include "Edges.hpp"
#include "MaxFlowDescription.hpp"

struct Net {
    std::vector<std::vector<InnerNetEdge> > graph;
    unsigned int source;
    unsigned int sink;
    unsigned int verticesCount;
    unsigned int edgesCount;
    
    Net();
    
    Net(unsigned int verticeCount, const std::vector<DirectedEdgeWithStart> &net, unsigned int source, unsigned int sink);
    
    MaxFlowDescription returnFlowDescription();
    
    void cleanUp();
};

#endif