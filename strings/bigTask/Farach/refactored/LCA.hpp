#ifndef _LCA
#define _LCA

#include <vector>
#include <algorithm>
#include <utility>
#include <climits>
#include <cstdio>
#include "sparseTable.hpp"
#include "RMQpm1.hpp"
#include "eulerPair.hpp"

using std::vector;
using std::max;
using std::min;

class LCA {
private:
    vector <unsigned int> myEuler;
    vector <unsigned int> first;
    vector <unsigned int> last;

    RMQpm1 rmq;

public:

    LCA();

    LCA(const vector<EulerPair> &euler);

    unsigned int lca(unsigned int u, unsigned int v) const;
};

#endif
