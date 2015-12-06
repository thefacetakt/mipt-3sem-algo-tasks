#ifndef _MINIMAL_PAIR
#define _MINIMAL_PAIR

template<typename T>
struct MinimalPair {
    T element;
    unsigned int index;

    MinimalPair() {
    }

    MinimalPair(T element, unsigned int index) :
        element(element), index(index) {
    }

    bool operator<(const MinimalPair &other) const {
        return (element == other.element ?
            index < other.index : element < other.element
        );
    }
};

#endif
