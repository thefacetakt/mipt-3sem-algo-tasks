#include "PrePushFlowSimpleFor.hpp"

#include "Edges.hpp"
#include "MaxFlowFinder.hpp"
#include <vector>
#include <algorithm>
#include <climits>

PrePushFlowSimpleFor::PrePushFlowSimpleFor() {
}

MaxFlowDescription PrePushFlowSimpleFor::findMaxFlow(unsigned int verticeCount, const std::vector<DirectedEdgeWithStart> &net, unsigned int source, unsigned int sink) {
    if (source == sink) {
        MaxFlowDescription result;
        result.flowValue = ULLONG_MAX;
        return result;
    }
    
    this->edgesCount = net.size();
    this->verticeCount = verticeCount;
    this->source = source;
    this->sink = sink;
    
    graph.resize(verticeCount);
    height.assign(verticeCount, 0);
    excess.assign(verticeCount, 0);
    firstUnsaturatedEdge.resize(verticeCount);
    
    for (unsigned int i = 0; i < net.size(); ++i) {
        const DirectedEdgeWithStart &e = net[i];
        
        graph[e.from].push_back(InnerNetEdge(e.to, e.capacity, graph[e.to].size(), i));
        graph[e.to].push_back(InnerNetEdge(e.from, 0, graph[e.from].size() - 1));
    }
    
    for (unsigned int v = 0; v < verticeCount; ++v) {
        firstUnsaturatedEdge[v] = graph[v].begin();
    }
    
    return findMaxFlowInitialised();
}

MaxFlowDescription PrePushFlowSimpleFor::findMaxFlowInitialised() {
    height[source] = verticeCount;
    
    unsigned int excessedCount = 0;

    for (auto &e: graph[source]) {
        e.flow = e.capacity;
        if (excess[e.to] == 0 && e.capacity)
            ++excessedCount;
        excess[e.to] += e.capacity;
        graph[e.to][e.reverse].flow -= e.capacity;
    }
    
    while(excessedCount) {
        for (unsigned int v = 0; v < verticeCount; ++v) {
            if (v != sink && v != source)
                discharge(v);
        }
        excessedCount = 0;
        for (unsigned int v = 0; v < verticeCount; ++v) {
            if (v != sink && v != source && excess[v])
                ++excessedCount;
        }
    }
    
    MaxFlowDescription result;
    result.flowValue = 0;
    for (auto const &e: graph[source]) {
        result.flowValue += e.flow;
    }
    
    result.description.resize(edgesCount);
    
    for (unsigned int v = 0; v < verticeCount; ++v) {
        for (auto const &e: graph[v]) {
            if (e.number != -1) {
                result.description[e.number] = FlowEdge(e.to, e.flow);
            }
        }
    }
    cleanUp();
    
    return result;
}

void PrePushFlowSimpleFor::discharge(unsigned int v) {
    while (excess[v]) {
        if (firstUnsaturatedEdge[v] == graph[v].end()) {
            relabel(v);
            firstUnsaturatedEdge[v] = graph[v].begin();
        } else {
            InnerNetEdge &e = *firstUnsaturatedEdge[v];
            if (e.capacity - e.flow > 0 && height[e.to] + 1 == height[v]) {
                push(v, e);
            } else {
                ++firstUnsaturatedEdge[v];
            }
        }
    }
}

void PrePushFlowSimpleFor::relabel(unsigned int v) {
    unsigned int currentHeight = UINT_MAX;
    for (auto const &e: graph[v]) {
        if (e.capacity - e.flow > 0) {
            currentHeight = std::min(currentHeight, height[e.to]);
        }
    }
    height[v] = currentHeight + 1;
}

void PrePushFlowSimpleFor::push(unsigned int v, InnerNetEdge &e) {
    long long pushValue = std::min(static_cast<unsigned long long> (e.capacity - e.flow), excess[v]);
    excess[v] -= pushValue;
    excess[e.to] += pushValue;
    e.flow += pushValue;
    graph[e.to][e.reverse].flow -= pushValue;
}

void PrePushFlowSimpleFor::cleanUp() {
    source = sink = verticeCount = edgesCount = 0;
    graph.clear();
    height.clear();
    excess.clear();
    firstUnsaturatedEdge.clear();
}