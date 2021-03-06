#ifndef _LCA_CPP
#define _LCA_CPP

#include <vector>
#include <algorithm>
#include <utility>
#include <climits>
#include <cstdio>
#include "sparseTable.hpp"
#include "RMQpm1.hpp"
#include "eulerPair.hpp"
#include "LCA.hpp"


using std::vector;
using std::max;
using std::min;

LCA::LCA() {
}

LCA::LCA(const vector<EulerPair> &euler) {
    unsigned int n = (euler.size() + 1) / 2;
    vector <unsigned int> toRMQ(euler.size());
    myEuler.resize(euler.size());
    first.assign(n, -1);
    last.assign(n, -1);
    for (unsigned int i = 0; i < euler.size(); ++i) {
        if (first[euler[i].node] == -1) {
            first[euler[i].node] = i;
        }
        last[euler[i].node] = i;
        toRMQ[i] = euler[i].depth;
        myEuler[i] = euler[i].node;
    }

    rmq = RMQpm1(toRMQ);
}


unsigned int LCA::lca(unsigned int u, unsigned int v) const {
    return myEuler[rmq.minimum(min(first[u], first[v]), max(last[u], last[v]))];
}

#endif
