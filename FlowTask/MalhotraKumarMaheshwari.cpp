#include <vector>
#include <queue>
#include "MalhotraKumarMaheshwari.hpp"
#include "Edges.hpp"
#include "MaxFlowFinder.hpp"
#include "MaxFlowDescription.hpp"

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

            
            for (auto const &e: graph[minimalNotBlocked]) {
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
        for (auto const &e: graph[v]) {
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
        for (auto const &e: graph[v]) {
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
