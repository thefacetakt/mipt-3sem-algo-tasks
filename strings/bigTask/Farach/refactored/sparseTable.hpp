#ifndef _SPARSE_TABLE
#define _SPARSE_TABLE


#include <vector>
#include <algorithm>
#include <climits>
#include "minimalPair.hpp"

using std::vector;
using std::min;

struct SparseTable {
private:
    vector<vector<MinimalPair> > st;
    vector<unsigned int> fastLog;
    unsigned int n;

    void init();

public:
    SparseTable(const vector<unsigned int> &elements);
    SparseTable();

    unsigned int minimum(unsigned int i, unsigned int j) const;

    MinimalPair operator[](unsigned int i) const;
};

#endif
