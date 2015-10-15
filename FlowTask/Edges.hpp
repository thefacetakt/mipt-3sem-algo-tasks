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
    
    InnerNetEdge(NetEdge edge, unsigned int reverse, unsigned int number=-1) : InnerNetEdge(edge.to, edge.capacity, reverse, number) {
    }
};


#endif
