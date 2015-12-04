#include <vector>
#include "Edges.hpp"
#include "Net.hpp"

Net::Net() {
}

Net::Net(unsigned int verticesCount, const std::vector<DirectedEdgeWithStart> &net, unsigned int source, unsigned int sink) {
    this->edgesCount = net.size();
    this->verticesCount = verticesCount;
    this->source = source;
    this->sink = sink;
    
    graph.resize(verticesCount);
    
    for (unsigned int i = 0; i < net.size(); ++i) {
        const DirectedEdgeWithStart &e = net[i];
        
        graph[e.from].push_back(InnerNetEdge(e.to, e.capacity, graph[e.to].size(), i));
        graph[e.to].push_back(InnerNetEdge(e.from, 0, graph[e.from].size() - 1));
    }
}
    
MaxFlowDescription Net::returnFlowDescription() {
    MaxFlowDescription result;
    result.flowValue = 0;
    for (auto const &e: graph[source]) {
        result.flowValue += e.flow;
    }
    
    result.description.resize(edgesCount);
    
    for (unsigned int v = 0; v < verticesCount; ++v) {
        for (auto const &e: graph[v]) {
            if (e.number != -1) {
                result.description[e.number] = FlowEdge(e.to, e.flow);
            }
        }
    }
    return result;
}

void Net::cleanUp() {
    source = sink = verticesCount = edgesCount = 0;
    graph.clear();
}