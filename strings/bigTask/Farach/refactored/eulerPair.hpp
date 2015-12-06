#ifndef _EULER_PAIR
#define _EULER_PAIR

struct EulerPair {
    unsigned int depth;
    unsigned int node;

    EulerPair();

    EulerPair(unsigned int depth, unsigned int node);
};

#endif
