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


RMQpm1::RMQpm1() {
}

RMQpm1::RMQpm1(const vector<unsigned int> &elements) {
    this->elements = elements;
    n = elements.size();
    for (block = 1; (1 << block) <= n; ++block) {
    }
    block = max(block / 2, 1u);
    vector <vector <unsigned int> > decomposition;
    vector <unsigned int> stElements;
    for (unsigned int i = 0; i < n; i += block) {
        decomposition.push_back(vector<unsigned int>(block, 0));
        prefixMins.push_back(vector <pair<unsigned int, unsigned int> >(block));
        suffixMins.push_back(vector <pair<unsigned int, unsigned int> >(block));
        
        unsigned int currentMask = 0;
        
        for (unsigned int j = 0; j < block; ++j) {
            if (i + j < n) {
                decomposition.back()[j] = elements[i + j];
            }
            if (j > 1 && decomposition.back()[j] > decomposition.back()[j - 1]) {
                currentMask |= (1 << j);
            }
        }
        currentMask >>= 1;
        type.push_back(currentMask);
        
        prefixMins.back()[0] = make_pair(decomposition.back()[0], i + 0);
        suffixMins.back()[block - 1] = make_pair(decomposition.back()[block - 1], i + block - 1);
        for (unsigned int j = 1; j < block; ++j) {
            prefixMins.back()[j] = min(prefixMins.back()[j - 1], make_pair(decomposition.back()[j], i + j));
            suffixMins.back()[block - j - 1] = min(suffixMins.back()[block - j], make_pair(decomposition.back()[block - j - 1], i + block - j - 1));
        }
        stElements.push_back(prefixMins.back().back().first);
    }
    st = SparseTable(stElements);

    dp.resize(1 << (block - 1));
    
    for (unsigned int i = 0; i < (1 << (block - 1)); ++i) {
        dp[i].resize(block);
        for (unsigned int length = 1; length <= block; ++length) {
            dp[i][length - 1].resize(block - length + 2);
            if (length == 1) {
                for (unsigned int j = 0; j < block; ++j) {
                    dp[i][length - 1][j] = make_pair((j == 0 ? 0 : dp[i][length - 1][j - 1].first + 2 * ((i >> (j - 1)) & 1) - 1), j);
                }
            } else {
                for (unsigned int j = length - 2; j < block; ++j) {
                    dp[i][length - 1][j - (length - 2)] = min(dp[i][0][j], dp[i][length - 2][j - (length - 2)]);
                }
            }
        }
    }
}
    
unsigned int RMQpm1::minimum(unsigned int i, unsigned int j) {
    if (j < i) {
        return UINT_MAX;
    }
    if (i == j) {
        return i;
    }
    unsigned int iBlock = i / block;
    unsigned int jBlock = j / block;
    if (iBlock != jBlock) {
        pair<unsigned int, unsigned int> prefSufMin = min(suffixMins[iBlock][i % block], prefixMins[jBlock][j % block]);
        if (iBlock + 1 > jBlock - 1) {
            return prefSufMin.second;
        } else {
            unsigned int insideMinPos = st.minimum(iBlock + 1, jBlock - 1);
            if (st[insideMinPos].first < prefSufMin.first) {
                return prefixMins[insideMinPos][block - 1].second;
            } else {
                return prefSufMin.second;
            }
        }
    } else {
        return iBlock * block + dp[type[iBlock]][j - i][i % block].second;
    }
}

#endif