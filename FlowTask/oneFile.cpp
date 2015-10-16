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

#ifndef _MALHOTA_KUMAR_MAHESHWARI_
#define _MALHOTA_KUMAR_MAHESHWARI_

#include <vector>
#include <climits>

class MalhotraKumarMaheshwari: public MaxFlowFinder {
    enum EDirections {
        FORWARD,
        BACKWARD,
        DIRECTIONS,
    };
    
    std::vector<std::vector<InnerNetEdge> > graph;
    std::vector<unsigned long long> potential[DIRECTIONS];
    std::vector<bool> blocked;
    std::vector<unsigned int> distance;
    std::vector<std::vector<InnerNetEdge>::iterator> firstUnsaturatedEdge[DIRECTIONS];
    
    
    unsigned int source;
    unsigned int sink;
    unsigned int verticesCount;
    unsigned int edgesCount;
    
    void cleanUp();
    
    unsigned int cf(const InnerNetEdge &e);
    
    InnerNetEdge &reverseEdge(const InnerNetEdge &e);
    
    void countAllPotentials();
    
    bool BFS();
    
    void pushVertice(unsigned int v, EDirections direction, unsigned long long value);
    
    void push(unsigned int v, InnerNetEdge &e, unsigned long long value);
    
    MaxFlowDescription findMaxFlowInitialised();
    
    unsigned long long phi(unsigned int v);
    
    
    
public:
    MalhotraKumarMaheshwari();
    
    MaxFlowDescription findMaxFlow(unsigned int verticeCount, const std::vector<DirectedEdgeWithStart> &net, unsigned int source, unsigned int sink);
    
    virtual ~MalhotraKumarMaheshwari() {
    }
    
};
#endif

#include <vector>
#include <queue>

MalhotraKumarMaheshwari::MalhotraKumarMaheshwari() {
}

MaxFlowDescription MalhotraKumarMaheshwari::findMaxFlow(unsigned int verticeCount, const std::vector<DirectedEdgeWithStart> &net, unsigned int source, unsigned int sink) {
    if (source == sink) {
        MaxFlowDescription result;
        result.flowValue = ULLONG_MAX;
        return result;
    }
    
    this->source = source;
    this->sink = sink;
    this->verticesCount = verticeCount;
    this->edgesCount = net.size();
    
    graph.resize(verticeCount);
    blocked.assign(verticeCount, false);
    for (unsigned int i = 0; i < net.size(); ++i) {
        const DirectedEdgeWithStart &e = net[i];
        
        graph[e.from].push_back(InnerNetEdge(e.to, e.capacity, graph[e.to].size(), i));
        graph[e.to].push_back(InnerNetEdge(e.from, 0, graph[e.from].size() - 1));
    }
        
    return findMaxFlowInitialised();
}

MaxFlowDescription MalhotraKumarMaheshwari::findMaxFlowInitialised() {
    while (BFS()) {
        countAllPotentials();
        bool notBlockedAvalible = true;
        while (notBlockedAvalible) {
            notBlockedAvalible = false;
            
            unsigned int minimalNotBlocked = 0;
            
            
            for (unsigned int v = 0; v < verticesCount; ++v) {
                if (phi(v) < phi(minimalNotBlocked)) {
                    minimalNotBlocked = v;
                }
            }

            unsigned long long pushValue = phi(minimalNotBlocked);
            
            pushVertice(minimalNotBlocked, FORWARD, pushValue);
            pushVertice(minimalNotBlocked, BACKWARD, pushValue);

            for (unsigned int i = 0; i < graph[minimalNotBlocked].size(); ++i) {
                const InnerNetEdge &e = graph[minimalNotBlocked][i];
                if (distance[e.to] == distance[minimalNotBlocked] + 1) {
                    potential[BACKWARD][e.to] -= cf(e);
                }
                if (distance[minimalNotBlocked] == distance[e.to] + 1) {
                    potential[FORWARD][e.to] -= cf(reverseEdge(e));
                }
            }
            
            blocked[minimalNotBlocked] = true;
            
            for (unsigned int v = 0; v < verticesCount; ++v) {
                if (!blocked[v]) {
                    notBlockedAvalible = true;
                }
            }
        }
    }
    MaxFlowDescription result;
    result.flowValue = 0;
    for (unsigned int i = 0; i < graph[source].size(); ++i) {
        const InnerNetEdge &e = graph[source][i];
        result.flowValue += e.flow;
    }
    
    result.description.resize(edgesCount);
    
    for (unsigned int v = 0; v < verticesCount; ++v) {
        for (unsigned int i = 0; i < graph[v].size(); ++i) {
            const InnerNetEdge &e = graph[v][i];
            if (e.number != -1) {
                result.description[e.number] = FlowEdge(e.to, e.flow);
            }
        }
    }
    cleanUp();
    return result;
}

void MalhotraKumarMaheshwari::countAllPotentials() {
    blocked.assign(verticesCount, false);
    for (unsigned int i = 0; i < DIRECTIONS; ++i) {
        potential[i].assign(verticesCount, 0);
        
        firstUnsaturatedEdge[i].resize(verticesCount);
        
        for (unsigned int v = 0; v < verticesCount; ++v) {
            firstUnsaturatedEdge[i][v] = graph[v].begin();
        }
    }
    
    for (unsigned int v = 0; v < verticesCount; ++v) {
        if (distance[v] >= distance[sink] && v != sink) {
            blocked[v] = true;
        }
    }
    
    for (unsigned int v = 0; v < verticesCount; ++v) {
        for (unsigned int i = 0; i < graph[v].size(); ++i) {
            const InnerNetEdge &e = graph[v][i];
            if (cf(e) > 0 && distance[e.to] == distance[v] + 1 && !blocked[e.to]) {
                potential[FORWARD][v] += cf(e);
                potential[BACKWARD][e.to] += cf(e);
            }
        }
    }
}

bool MalhotraKumarMaheshwari::BFS() {
    std::queue<unsigned int> q;
    q.push(source);
    distance.assign(verticesCount, UINT_MAX);
    distance[source] = 0;
    while (!q.empty()) {
        unsigned int v = q.front();
        q.pop();
        for (unsigned int i = 0; i < graph[v].size(); ++i) {
            const InnerNetEdge &e = graph[v][i];
            if (cf(e) > 0 && distance[e.to] > distance[v] + 1) {
                distance[e.to] = distance[v] + 1;
                q.push(e.to);
            }
        }
    }
    return (distance[sink] != UINT_MAX);
}

unsigned long long MalhotraKumarMaheshwari::phi(unsigned int v) {
    if (blocked[v])
        return ULLONG_MAX;
    if (v == source)
        return potential[FORWARD][v];
    if (v == sink)
        return potential[BACKWARD][v];
    
    return std::min(potential[FORWARD][v], potential[BACKWARD][v]);
}

void MalhotraKumarMaheshwari::pushVertice(unsigned int v, EDirections direction, unsigned long long value) {
    std::vector<unsigned long long> pushValue(verticesCount, 0);
    std::queue<unsigned int> touched;
    touched.push(v);
    pushValue[v] = value;
    while (!touched.empty()) {
        unsigned int currentV = touched.front();
        touched.pop();
        for (; pushValue[currentV] && firstUnsaturatedEdge[direction][currentV] != graph[currentV].end(); ++firstUnsaturatedEdge[direction][currentV]) {
            InnerNetEdge &e = *firstUnsaturatedEdge[direction][currentV];
            if (distance[e.to] == distance[currentV] + 1 - 2 * direction && !blocked[e.to]) {
                unsigned int u = currentV;
                InnerNetEdge *candidateEdge = &e;
                if (direction == BACKWARD) {
                    u = e.to;
                    candidateEdge = &reverseEdge(e);
                }
                InnerNetEdge &pushEdge = *candidateEdge;
                if (cf(pushEdge) > 0) {
                    unsigned long long value = std::min(static_cast<unsigned long long>(cf(pushEdge)), pushValue[currentV]);
                    push(u, pushEdge, value);
                    pushValue[currentV] -= value;
                    pushValue[e.to] += value;
                    touched.push(e.to);
                    if (!pushValue[currentV])
                        break;
                }
            }
        }
    }
}

void MalhotraKumarMaheshwari::push(unsigned int v, InnerNetEdge &e, unsigned long long value) {
    e.flow += value;
    reverseEdge(e).flow -= value;
    potential[FORWARD][v] -= value;
    potential[BACKWARD][e.to] -= value;
}

void MalhotraKumarMaheshwari::cleanUp() {
    source = sink = verticesCount = edgesCount = 0;
    graph.clear();
    blocked.clear();
    distance.clear();
    for (unsigned int i = 0; i < DIRECTIONS; ++i) {
        potential[i].clear();
        firstUnsaturatedEdge[i].clear();
    }
}

unsigned int MalhotraKumarMaheshwari::cf(const InnerNetEdge &e) {
    return e.capacity - e.flow;
}

InnerNetEdge & MalhotraKumarMaheshwari::reverseEdge(const InnerNetEdge &e) {
    return graph[e.to][e.reverse];
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
                return new MalhotraKumarMaheshwari();
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
    MaxFlowFinder *maxFlowFinder = MaxFlowFinderFabric::getMaxFlowFinder(MaxFlowFinderFabric::MALHOTRA_KUMAR_MAHESHWARI);
    MaxFlowDescription description = maxFlowFinder->findMaxFlow(n, graph, 0, n - 1);
    printf("%llu\n",  description.flowValue);
    for (unsigned int i = 0; i < description.description.size(); ++i) {
        printf("%u\n", description.description[i].flow);
    }
    delete maxFlowFinder;
    return 0;
}