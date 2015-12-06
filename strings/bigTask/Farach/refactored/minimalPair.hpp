#ifndef _MINIMAL_PAIR
#define _MINIMAL_PAIR

struct MinimalPair {
    unsigned int element;
    unsigned int index;

    MinimalPair();

    MinimalPair(unsigned int element, unsigned int index);

    bool operator<(const MinimalPair &other) const;
};

#endif
