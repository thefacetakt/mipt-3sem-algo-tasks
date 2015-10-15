#ifndef _EDGE_
#define _EDGE_

class DirectedEdgeWithStart {
public:
    unsigned int from;
    unsigned int to;
    unsigned int capacity;
    
    DirectedEdgeWithStart() {
    }
    
    DirectedEdgeWithStart(unsigned int from, unsigned int to, unsigned int capacity) : from(from), to(to), capacity(capacity) {
    }
};


class Edge {
public:
    unsigned int to;
    
    Edge() {
    }
    
    Edge(unsigned int to) : to(to) {
    }
};

class NetEdge: public Edge {
public:
    unsigned int capacity;
    NetEdge(unsigned int to, unsigned int capacity) : Edge(to), capacity(capacity) {
    }
};

class FlowEdge: public Edge {
public:
    int flow;
    
    FlowEdge() {
    }
    
    FlowEdge(unsigned int to, int flow) : Edge(to), flow(flow) {
    }
};

class InnerNetEdge: public NetEdge {
public:
    int flow;
    unsigned int reverse;
    unsigned int number;
    
    InnerNetEdge(unsigned int to, unsigned int capacity, unsigned int reverse, unsigned int number=-1) : NetEdge(to, capacity), flow(0), reverse(reverse), number(number) {
    }
};


#endif



#ifndef _MAX_FLOW_DESCRIPTION_
#define _MAX_FLOW_DESCRIPTION_


#include <vector>

class MaxFlowDescription {
public:
    unsigned long long flowValue;
    std::vector<FlowEdge> description;
};

#endif


#ifndef _MAX_FLOW_FINDER_
#define _MAX_FLOW_FINDER_


#include <vector>

class MaxFlowFinder {
public:
    virtual MaxFlowDescription findMaxFlow(unsigned int verticeCount, const std::vector<DirectedEdgeWithStart> &net, unsigned int source, unsigned int sink) = 0;
    virtual ~MaxFlowFinder() {
    }
};



#endif

#ifndef _PRE_PUSH_FLOW_SIMPLE_FOR_
#define _PRE_PUSH_FLOW_SIMPLE_FOR_

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

#include <vector>
#include <algorithm>
#include <climits>

PrePushFlowSimpleFor::PrePushFlowSimpleFor() {
}

MaxFlowDescription PrePushFlowSimpleFor::findMaxFlow(unsigned int verticeCount, const std::vector<DirectedEdgeWithStart> &net, unsigned int source, unsigned int sink) {
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

    for (int i = 0; i < graph[source].size(); ++i) {
        InnerNetEdge &e = graph[source][i];
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
    for (int i = 0; i < graph[source].size(); ++i) {
        const InnerNetEdge &e = graph[source][i];
        result.flowValue += e.flow;
    }
    
    result.description.resize(edgesCount);
    
    for (unsigned int v = 0; v < verticeCount; ++v) {
        for (int i = 0; i < graph[v].size(); ++i) {
        const InnerNetEdge &e = graph[v][i];
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
    for (int i = 0; i < graph[v].size(); ++i) {
        const InnerNetEdge &e = graph[v][i];
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

#ifndef _MAX_FLOW_FINDER_FABRIC_
#define _MAX_FLOW_FINDER_FABRIC_


class MaxFlowFinderFabric {
public:
    enum EAlgorithms {
        PRE_PUSH_FLOW_SIMPLE_FOR,
        MALHOTRA_KUMAR_MAHESHWARI,
    };
    
    static MaxFlowFinder *getMaxFlowFinder(EAlgorithms algorithm) {
        switch(algorithm) {
            case PRE_PUSH_FLOW_SIMPLE_FOR:
                return new PrePushFlowSimpleFor();
            case MALHOTRA_KUMAR_MAHESHWARI:
                return NULL; //new MalhotraKumarMaheshwari();
        }
    };
};

#endif

#include <cstdio>
#include <vector>

int main() {
    unsigned int n, m;
    scanf("%u %u", &n, &m);
    std::vector <DirectedEdgeWithStart> graph(m);
    for (unsigned int i = 0; i < m; ++i) {
        unsigned int from, to, capacity;
        scanf("%u %u %u", &from, &to, &capacity);
        --from, --to;
        graph[i] = DirectedEdgeWithStart(from, to, capacity);
    }
    MaxFlowFinder *maxFlowFinder = MaxFlowFinderFabric::getMaxFlowFinder(MaxFlowFinderFabric::PRE_PUSH_FLOW_SIMPLE_FOR);
    MaxFlowDescription description = maxFlowFinder->findMaxFlow(n, graph, 0, n - 1);
    printf("%llu\n",  description.flowValue);
    for (unsigned int i = 0; i < description.description.size(); ++i) {
        printf("%u\n", description.description[i].flow);
    }
    delete maxFlowFinder;
    return 0;
}