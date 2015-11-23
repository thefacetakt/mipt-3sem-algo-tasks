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
    
    LCA(const vector<pair<unsigned int, unsigned int> > &euler);
    
    
    unsigned int lca(unsigned int u, unsigned int v);
};
    
LCA::LCA(const vector<pair<unsigned int, unsigned int> > &euler) {
    unsigned int n = (euler.size() + 1) / 2;
    vector <unsigned int> toRMQ(euler.size());
    myEuler.resize(euler.size());
    first.assign(n, -1);
    last.assign(n, -1);
    for (unsigned int i = 0; i < euler.size(); ++i) {
        if (first[euler[i].second] == -1) {
            first[euler[i].second] = i;
        }
        last[euler[i].second] = i;
        toRMQ[i] = euler[i].first;
        myEuler[i] = euler[i].second;
    }
    
    
    rmq = RMQpm1(toRMQ);
}


unsigned int LCA::lca(unsigned int u, unsigned int v) {
    return myEuler[rmq.minimum(min(first[u], first[v]), max(last[u], last[v]))];
}


#endif