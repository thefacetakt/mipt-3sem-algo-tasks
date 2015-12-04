#include "PrePushFlowSimpleFor.hpp"

#include "Edges.hpp"
#include "MaxFlowFinder.hpp"
#include <vector>
#include <algorithm>
#include <climits>

PrePushFlowSimpleFor::PrePushFlowSimpleFor() {
}

MaxFlowDescription PrePushFlowSimpleFor::findMaxFlow(unsigned int verticesCount, const std::vector<DirectedEdgeWithStart> &net, unsigned int source, unsigned int sink) {
    if (source == sink) {
        MaxFlowDescription result;
        result.flowValue = ULLONG_MAX;
        return result;
    }
    
    net_ = Net(verticesCount, net, source, sink);
    
    height_.assign(verticesCount, 0);
    excess_.assign(verticesCount, 0);
    firstUnsaturatedEdge_.resize(verticesCount);
    
    for (unsigned int v = 0; v < verticesCount; ++v) {
        firstUnsaturatedEdge_[v] = net_.graph[v].begin();
    }
    
    return findMaxFlowInitialised();
}

MaxFlowDescription PrePushFlowSimpleFor::findMaxFlowInitialised() {
    height_[net_.source] = net_.verticesCount;
    
    unsigned int excessedCount = 0;

    for (auto &e: net_.graph[net_.source]) {
        e.flow = e.capacity;
        if (excess_[e.to] == 0 && e.capacity)
            ++excessedCount;
        excess_[e.to] += e.capacity;
        net_.graph[e.to][e.reverse].flow -= e.capacity;
    }
    
    while(excessedCount) {
        for (unsigned int v = 0; v < net_.verticesCount; ++v) {
            if (v != net_.sink && v != net_.source)
                discharge(v);
        }
        excessedCount = 0;
        for (unsigned int v = 0; v < net_.verticesCount; ++v) {
            if (v != net_.sink && v != net_.source && excess_[v])
                ++excessedCount;
        }
    }
    
    MaxFlowDescription result = net_.returnFlowDescription();
    cleanUp();
    
    return result;
}

void PrePushFlowSimpleFor::discharge(unsigned int v) {
    while (excess_[v]) {
        if (firstUnsaturatedEdge_[v] == net_.graph[v].end()) {
            relabel(v);
            firstUnsaturatedEdge_[v] = net_.graph[v].begin();
        } else {
            InnerNetEdge &e = *firstUnsaturatedEdge_[v];
            if (e.capacity - e.flow > 0 && height_[e.to] + 1 == height_[v]) {
                push(v, e);
            } else {
                ++firstUnsaturatedEdge_[v];
            }
        }
    }
}

void PrePushFlowSimpleFor::relabel(unsigned int v) {
    unsigned int currentHeight = UINT_MAX;
    for (auto const &e: net_.graph[v]) {
        if (e.capacity - e.flow > 0) {
            currentHeight = std::min(currentHeight, height_[e.to]);
        }
    }
    height_[v] = currentHeight + 1;
}

void PrePushFlowSimpleFor::push(unsigned int v, InnerNetEdge &e) {
    long long pushValue = std::min(static_cast<unsigned long long> (e.capacity - e.flow), excess_[v]);
    excess_[v] -= pushValue;
    excess_[e.to] += pushValue;
    e.flow += pushValue;
    net_.graph[e.to][e.reverse].flow -= pushValue;
}

void PrePushFlowSimpleFor::cleanUp() {
    net_.cleanUp();
    height_.clear();
    excess_.clear();
    firstUnsaturatedEdge_.clear();
}