#include <vector>
#include <algorithm>
#include <climits>
#include <cstdio>
#include "sparseTable.hpp"
#include "RMQpm1.hpp"
#include "minimalPair.hpp"

using std::vector;
using std::max;
using std::min;

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
        prefixMins.push_back(vector <MinimalPair<unsigned int> >(block));
        suffixMins.push_back(vector <MinimalPair<unsigned int> >(block));

        unsigned int currentMask = 0;

        for (unsigned int j = 0; j < block; ++j) {
            if (i + j < n) {
                decomposition.back()[j] = elements[i + j];
            }
            if (j >= 1
                && decomposition.back()[j] > decomposition.back()[j - 1]) {
                currentMask |= (1 << j);
            }
        }
        currentMask >>= 1;
        type.push_back(currentMask);

        prefixMins.back()[0]
            = MinimalPair<unsigned int>(decomposition.back()[0], i + 0);
        suffixMins.back()[block - 1]
            = MinimalPair<unsigned int>(decomposition.back()[block - 1],
                i + block - 1);
        for (unsigned int j = 1; j < block; ++j) {
            prefixMins.back()[j] = min(prefixMins.back()[j - 1],
                MinimalPair<unsigned int>(decomposition.back()[j], i + j)
            );
            suffixMins.back()[block - j - 1] = min(suffixMins.back()[block - j],
                MinimalPair<unsigned int>(decomposition.back()[block - j - 1],
                    i + block - j - 1
                )
            );
        }
        stElements.push_back(prefixMins.back().back().element);
    }
    st = SparseTable(stElements);

    dp.resize(1 << (block - 1));

    for (unsigned int i = 0; i < (1 << (block - 1)); ++i) {
        dp[i].resize(block);
        for (unsigned int length = 1; length <= block; ++length) {
            dp[i][length - 1].resize(block - length + 2);
            if (length == 1) {
                for (unsigned int j = 0; j < block; ++j) {
                    dp[i][length - 1][j] = MinimalPair<int>((j == 0 ?
                        0 : dp[i][length - 1][j - 1].element
                            + 2 * ((i >> (j - 1)) & 1) - 1
                        ), j
                    );
                }
            } else {
                for (unsigned int j = length - 2; j < block; ++j) {
                    dp[i][length - 1][j - (length - 2)] = min(dp[i][0][j],
                        dp[i][length - 2][j - (length - 2)]
                    );
                }
            }
        }
    }
}

unsigned int RMQpm1::minimum(unsigned int i, unsigned int j) const {
    if (j < i) {
        return UINT_MAX;
    }
    if (i == j) {
        return i;
    }
    unsigned int iBlock = i / block;
    unsigned int jBlock = j / block;
    if (iBlock != jBlock) {
        MinimalPair<unsigned int> prefSufMin
            = min(suffixMins[iBlock][i % block], prefixMins[jBlock][j % block]);
        if (iBlock + 1 > jBlock - 1) {
            return prefSufMin.index;
        } else {
            unsigned int insideMinPos = st.minimum(iBlock + 1, jBlock - 1);
            if (st[insideMinPos].element < prefSufMin.element) {
                return prefixMins[insideMinPos][block - 1].index;
            } else {
                return prefSufMin.index;
            }
        }
    } else {
        return iBlock * block + dp[type[iBlock]][j - i][i % block].index;
    }
}
