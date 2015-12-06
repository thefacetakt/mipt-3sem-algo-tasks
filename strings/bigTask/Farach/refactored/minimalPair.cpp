#include "minimalPair.hpp"

MinimalPair::MinimalPair() {
}

MinimalPair::MinimalPair(unsigned int element, unsigned int index) :
    element(element),index(index) {
}

bool MinimalPair::operator<(const MinimalPair &other) const {
    return (element == other.element ?
        index < other.index : element < other.element
    );
}
