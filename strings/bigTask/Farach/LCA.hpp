#ifndef _LCA
#define _LCA

#include <vector>
#include <algorithm>
#include <utility>
#include <climits>
#include <cstdio>
#include "sparseTable.hpp"
#include "RMQpm1.hpp"

using std::vector;
using std::max;
using std::min;
using std::pair;
using std::make_pair;

struct LCA {
private:
    vector <unsigned int> myEuler;
    vector <unsigned int> first;
    vector <unsigned int> last;
    
    RMQpm1 rmq;
    
public:
    
    LCA();
    
    LCA(const vector<pair<unsigned int, unsigned int> > &euler);
    
    
    unsigned int lca(unsigned int u, unsigned int v);
};

#endif