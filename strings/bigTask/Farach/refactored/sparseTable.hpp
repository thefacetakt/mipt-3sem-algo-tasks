#ifndef _SPARSE_TABLE
#define _SPARSE_TABLE


#include <vector>
#include <algorithm>
#include <climits>
#include <utility>

using std::vector;
using std::min;
using std::pair;
using std::make_pair;

struct SparseTable {
private:
    vector<vector<pair<unsigned int, unsigned int> > > st;
    vector<unsigned int> fastLog;
    unsigned int n;

    void init();

public:
    SparseTable(const vector<unsigned int> &elements);
    SparseTable();

    unsigned int minimum(unsigned int i, unsigned int j) const;

    pair<unsigned int, unsigned int> operator[](unsigned int i) const {
        return st[0][i];
    }
};

#endif
