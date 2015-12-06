#ifndef _USEFUL_STRUCTURES_IS
#define _USEFUL_STRUCTURES_IS

#include <vector>

using std::vector;

enum Type {
    MINUS,
    PLUS,
    STAR,
    TOTAL
};

struct Node {
    vector <unsigned int> children;
    vector <unsigned int> term;
};

struct StringSortingItem {
    unsigned int position;
    unsigned int string;

    StringSortingItem();

    StringSortingItem(unsigned int position, unsigned int string);
};

#endif
