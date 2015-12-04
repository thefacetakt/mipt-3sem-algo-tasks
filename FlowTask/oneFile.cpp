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

#ifndef _NET
#define _NET

#include <vector>

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

#include <vector>

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
    for (unsigned int z = 0; z < graph[source].size(); ++z) {
        InnerNetEdge &e = graph[source][z];
        result.flowValue += e.flow;
    }
    
    result.description.resize(edgesCount);
    
    for (unsigned int v = 0; v < verticesCount; ++v) {
        for (unsigned int z = 0; z < graph[v].size(); ++z) {
            InnerNetEdge &e = graph[v][z];
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




#ifndef _PRE_PUSH_FLOW_SIMPLE_FOR_
#define _PRE_PUSH_FLOW_SIMPLE_FOR_

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
    for (unsigned int z = 0; z < net_.graph[net_.source].size(); ++z) {
        InnerNetEdge &e = net_.graph[net_.source][z];
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
    for (unsigned int z = 0; z < net_.graph[v].size(); ++z) {
        InnerNetEdge &e = net_.graph[v][z];
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
    
    Net net_;
    
    
    std::vector<unsigned long long> potential_[DIRECTIONS];
    std::vector<bool> blocked_;
    std::vector<unsigned int> distance_;
    std::vector<std::vector<InnerNetEdge>::iterator> firstUnsaturatedEdge_[DIRECTIONS];
    
    void cleanUp();
    
    unsigned int residue(const InnerNetEdge &e);
    
    InnerNetEdge &reverseEdge(const InnerNetEdge &e);
    
    void countAllPotentials();
    
    bool BFS();
    
    void pushVertice(unsigned int v, EDirections direction, unsigned long long value);
    
    void push(unsigned int v, InnerNetEdge &e, unsigned long long value);
    
    MaxFlowDescription findMaxFlowInitialised();
    
    unsigned long long phi(unsigned int v);
    
    
    
public:
    MalhotraKumarMaheshwari();
    
    MaxFlowDescription findMaxFlow(unsigned int verticesCount, const std::vector<DirectedEdgeWithStart> &net, unsigned int source, unsigned int sink);
    
    virtual ~MalhotraKumarMaheshwari() {
    }
    
};
#endif

#include <vector>
#include <queue>

MalhotraKumarMaheshwari::MalhotraKumarMaheshwari() {
}

MaxFlowDescription MalhotraKumarMaheshwari::findMaxFlow(unsigned int verticesCount, const std::vector<DirectedEdgeWithStart> &net, unsigned int source, unsigned int sink) {
    if (source == sink) {
        MaxFlowDescription result;
        result.flowValue = ULLONG_MAX;
        return result;
    }
    net_ = Net(verticesCount, net, source, sink);
    
    blocked_.assign(verticesCount, false);
        
    return findMaxFlowInitialised();
}

MaxFlowDescription MalhotraKumarMaheshwari::findMaxFlowInitialised() {
    while (BFS()) {
        countAllPotentials();
        bool notBlockedAvalible = true;
        while (notBlockedAvalible) {
            notBlockedAvalible = false;
            
            unsigned int minimalNotBlocked = 0;
            
            
            for (unsigned int v = 0; v < net_.verticesCount; ++v) {
                if (phi(v) < phi(minimalNotBlocked)) {
                    minimalNotBlocked = v;
                }
            }

            unsigned long long pushValue = phi(minimalNotBlocked);
            
            pushVertice(minimalNotBlocked, FORWARD, pushValue);
            pushVertice(minimalNotBlocked, BACKWARD, pushValue);
    
            for (unsigned int z = 0; z < net_.graph[minimalNotBlocked].size(); ++z) {
                InnerNetEdge &e = net_.graph[minimalNotBlocked][z];
                if (distance_[e.to] == distance_[minimalNotBlocked] + 1) {
                    potential_[BACKWARD][e.to] -= residue(e);
                }
                if (distance_[minimalNotBlocked] == distance_[e.to] + 1) {
                    potential_[FORWARD][e.to] -= residue(reverseEdge(e));
                }
            }
            
            blocked_[minimalNotBlocked] = true;
            
            for (unsigned int v = 0; v < net_.verticesCount; ++v) {
                if (!blocked_[v]) {
                    notBlockedAvalible = true;
                }
            }
        }
    }
    MaxFlowDescription result = net_.returnFlowDescription();
    
    cleanUp();
    
    return result;
}

void MalhotraKumarMaheshwari::countAllPotentials() {
    blocked_.assign(net_.verticesCount, false);
    for (unsigned int i = 0; i < DIRECTIONS; ++i) {
        potential_[i].assign(net_.verticesCount, 0);
        
        firstUnsaturatedEdge_[i].resize(net_.verticesCount);
        
        for (unsigned int v = 0; v < net_.verticesCount; ++v) {
            firstUnsaturatedEdge_[i][v] = net_.graph[v].begin();
        }
    }
    
    for (unsigned int v = 0; v < net_.verticesCount; ++v) {
        if (distance_[v] >= distance_[net_.sink] && v != net_.sink) {
            blocked_[v] = true;
        }
    }
    
    for (unsigned int v = 0; v < net_.verticesCount; ++v) {
        for (unsigned int z = 0; z < net_.graph[v].size(); ++z) {
            InnerNetEdge &e = net_.graph[v][z];
            if (residue(e) > 0 && distance_[e.to] == distance_[v] + 1 && !blocked_[e.to]) {
                potential_[FORWARD][v] += residue(e);
                potential_[BACKWARD][e.to] += residue(e);
            }
        }
    }
}

bool MalhotraKumarMaheshwari::BFS() {
    std::queue<unsigned int> q;
    q.push(net_.source);
    distance_.assign(net_.verticesCount, UINT_MAX);
    distance_[net_.source] = 0;
    while (!q.empty()) {
        unsigned int v = q.front();
        q.pop();
        for (unsigned int z = 0; z < net_.graph[v].size(); ++z) {
            InnerNetEdge &e = net_.graph[v][z];
            if (residue(e) > 0 && distance_[e.to] > distance_[v] + 1) {
                distance_[e.to] = distance_[v] + 1;
                q.push(e.to);
            }
        }
    }
    return (distance_[net_.sink] != UINT_MAX);
}

unsigned long long MalhotraKumarMaheshwari::phi(unsigned int v) {
    if (blocked_[v])
        return ULLONG_MAX;
    if (v == net_.source)
        return potential_[FORWARD][v];
    if (v == net_.sink)
        return potential_[BACKWARD][v];
    
    return std::min(potential_[FORWARD][v], potential_[BACKWARD][v]);
}

void MalhotraKumarMaheshwari::pushVertice(unsigned int v, EDirections direction, unsigned long long value) {
    std::vector<unsigned long long> pushValue(net_.verticesCount, 0);
    std::queue<unsigned int> touched;
    touched.push(v);
    pushValue[v] = value;
    while (!touched.empty()) {
        unsigned int currentV = touched.front();
        touched.pop();
        for (; pushValue[currentV] && firstUnsaturatedEdge_[direction][currentV] != net_.graph[currentV].end(); ++firstUnsaturatedEdge_[direction][currentV]) {
            InnerNetEdge &e = *firstUnsaturatedEdge_[direction][currentV];
            if (distance_[e.to] == distance_[currentV] + 1 - 2 * direction && !blocked_[e.to]) {
                unsigned int u = currentV;
                InnerNetEdge *candidateEdge = &e;
                if (direction == BACKWARD) {
                    u = e.to;
                    candidateEdge = &reverseEdge(e);
                }
                InnerNetEdge &pushEdge = *candidateEdge;
                if (residue(pushEdge) > 0) {
                    unsigned long long value = std::min(static_cast<unsigned long long>(residue(pushEdge)), pushValue[currentV]);
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
    potential_[FORWARD][v] -= value;
    potential_[BACKWARD][e.to] -= value;
}

void MalhotraKumarMaheshwari::cleanUp() {
    net_.cleanUp();
    blocked_.clear();
    distance_.clear();
    for (unsigned int i = 0; i < DIRECTIONS; ++i) {
        potential_[i].clear();
        firstUnsaturatedEdge_[i].clear();
    }
}

unsigned int MalhotraKumarMaheshwari::residue(const InnerNetEdge &e) {
    return e.capacity - e.flow;
}

InnerNetEdge & MalhotraKumarMaheshwari::reverseEdge(const InnerNetEdge &e) {
    return net_.graph[e.to][e.reverse];
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
    MaxFlowFinder *maxFlowFinder = MaxFlowFinderFabric::getMaxFlowFinder(MaxFlowFinderFabric::PRE_PUSH_FLOW_SIMPLE_FOR);
    MaxFlowDescription description = maxFlowFinder->findMaxFlow(n, graph, 0, n - 1);
    printf("%llu\n",  description.flowValue);
    for (unsigned int i = 0; i < description.description.size(); ++i) {
        printf("%u\n", description.description[i].flow);
    }
    delete maxFlowFinder;
    return 0;
}