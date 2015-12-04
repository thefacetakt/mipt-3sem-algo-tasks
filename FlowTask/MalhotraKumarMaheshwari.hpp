#ifndef _MALHOTRA_KUMAR_MAHESHWARI_
#define _MALHOTRA_KUMAR_MAHESHWARI_

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
    
    unsigned long long phi(unsigned int v);
    
    void init(unsigned int verticesCount);
    
public:
    MalhotraKumarMaheshwari();
    
    void findMaxFlow();
    
    virtual ~MalhotraKumarMaheshwari() {
    }
    
};
#endif