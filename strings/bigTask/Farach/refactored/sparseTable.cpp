#include "sparseTable.hpp"

#include <vector>
#include <algorithm>
#include <climits>
#include <utility>

using std::vector;
using std::min;
using std::pair;
using std::make_pair;


SparseTable::SparseTable() {
}

SparseTable::SparseTable(const vector<unsigned int> &elements) {
    n = elements.size();
    fastLog.resize(n + 1);
    fastLog[1] = 0;
    for (unsigned int i = 2; i <= n; ++i) {
        fastLog[i] = fastLog[i / 2] + 1;
    }
    st.resize(fastLog[n] + 1);
    st[0].resize(n);
    for (unsigned int i = 0; i < n; ++i) {
        st[0][i] = make_pair(elements[i], i);
    }
    init();
}

void SparseTable::init() {
    for (unsigned int k = 1; k <= fastLog[n]; ++k) {
        st[k].resize(n - (1 << k) + 1);
        for (unsigned int i = 0; i < st[k].size(); ++i) {
            st[k][i] = min(st[k - 1][i], st[k - 1][i + (1 << (k - 1))]);
        }
    }
}

unsigned int SparseTable::minimum(unsigned int i, unsigned int j) const {
    if (j < i) {
        return UINT_MAX;
    }
    unsigned int length = (j - i + 1);
    return min(st[fastLog[length]][i], st[fastLog[length]][j - (1 << fastLog[length]) + 1]).second;
}
