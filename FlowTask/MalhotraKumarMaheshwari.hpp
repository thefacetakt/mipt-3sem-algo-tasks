#ifndef _MALHOTA_KUMAR_MAHESHWARI_
#define _MALHOTA_KUMAR_MAHESHWARI_

#include <vector>
#include <climits>
#include "Edges.hpp"
#include "MaxFlowFinder.hpp"
#include "MaxFlowDescription.hpp"

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