#ifndef _RMQpm1
#define _RMQpm1

#include <vector>
#include <algorithm>
#include <climits>
#include <cstdio>
#include "sparseTable.hpp"
#include "minimalPair.hpp"

using std::vector;
using std::max;
using std::min;

struct RMQpm1 {
private:
    SparseTable st;
    vector <unsigned int> elements;
    vector <vector <MinimalPair> > prefixMins;
    vector <vector <MinimalPair> > suffixMins;
    unsigned int n;
    unsigned int block;
    vector <vector<vector<MinimalPair> > > dp;
    vector <unsigned int> type;

public:
    RMQpm1();
    RMQpm1(const vector<unsigned int> &elements);

    unsigned int minimum(unsigned int i, unsigned int j) const;
};

#endif
