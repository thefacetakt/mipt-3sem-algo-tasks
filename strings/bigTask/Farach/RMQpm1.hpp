#ifndef _RMQpm1
#define _RMQpm1

#include <vector>
#include <algorithm>
#include <utility>
#include <climits>
#include <cstdio>
#include "sparseTable.hpp"

using std::vector;
using std::max;
using std::min;
using std::pair;
using std::make_pair;

struct RMQpm1 {
private:
    SparseTable st;
    vector <unsigned int> elements;
    vector <vector <pair<unsigned int, unsigned int> > > prefixMins;
    vector <vector <pair<unsigned int, unsigned int> > > suffixMins;
    unsigned int n;
    unsigned int block;
    vector <vector<vector<pair<int, unsigned int> > > > dp;
    vector <unsigned int> type;
    
public:
    RMQpm1();
    RMQpm1(const vector<unsigned int> &elements);
    
    unsigned int minimum(unsigned int i, unsigned int j);
};

#endif