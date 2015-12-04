#include <vector>
#include <queue>
#include "Net.hpp"
#include "MalhotraKumarMaheshwari.hpp"
#include "Edges.hpp"
#include "MaxFlowFinder.hpp"
#include "MaxFlowDescription.hpp"

MalhotraKumarMaheshwari::MalhotraKumarMaheshwari() {
}

void MalhotraKumarMaheshwari::init(unsigned int verticesCount) {
    blocked_.assign(verticesCount, false);
}

void MalhotraKumarMaheshwari::findMaxFlow() {
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

            
            for (auto const &e: net_.graph[minimalNotBlocked]) {
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
        for (auto const &e: net_.graph[v]) {
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
        for (auto const &e: net_.graph[v]) {
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
