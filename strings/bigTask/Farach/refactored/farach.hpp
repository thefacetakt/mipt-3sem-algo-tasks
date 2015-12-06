#ifndef _FARACH
#define _FARACH

#define _GLIBCXX_DEBUG
#define _PRINT_DBG
#define _DEBUG

#include <cstdio>
#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>
#include <climits>
#include <cassert>
#include <functional>
#include <string>
#include "EulerPair.hpp"
#include "LCA.hpp"
#include "SuffixTree.hpp"
#include "usefulStructures.hpp"

using std::vector;
using std::make_pair;
using std::max;
using std::min;
using std::pair;
using std::function;
using std::swap;

void buildEulerDfs(unsigned int v, unsigned int &count, unsigned int depth,
    const SuffixTree &tree, vector <EulerPair> &euler,
    vector <unsigned int> &realNode
);

void buildEuler(const SuffixTree &tree, vector <EulerPair> &euler,
    vector <unsigned int> &realNode
);



RandomLCPGetter usualGetter(const SuffixTree &tree, unsigned int inputLength);


SuffixTree buildTempSuffixTree(const vector <int> &);

template<class T, class Compare>
vector <T> radixSort(const vector <T> &input, Compare comp) {
    int maxElement = 0;
    for (auto const &element: input) {
        maxElement = max(maxElement, static_cast<int>(comp(element) + 1));
    }
    vector <unsigned int> count(maxElement + 1);
    for (auto const &element: input) {
        ++count[comp(element) + 1];
    }
    for (unsigned int i = 1; i <= maxElement; ++i) {
        count[i] += count[i - 1];
    }
    vector <T> result(input.size());
    for (unsigned int i = input.size() - 1; i != UINT_MAX; --i) {
        const T &element = input[i];
        result[--count[comp(element) + 1]] = element;
    }
    return result;
}

int firstElementOfNumberedPair(const NumberedPair &x);

int secondElementOfNumberedPair(const NumberedPair &x);

vector <int> compressInput(const vector <int> &input);

int decompressDfs(unsigned int v, unsigned int inParentIndex, SuffixTree &tree,
    const vector <int> &input, unsigned int depth
);

void decompress(SuffixTree &tree, const vector <int> &input);

void checkDfs(const SuffixTree &tree, unsigned int v, const vector<int> &input);

void printHiddenDfs(const SuffixTree &tree, unsigned int v);

void suffixArrayDfs(unsigned int v, const SuffixTree &tree,
    vector <unsigned int> &suffix);

void buildSuffixArray(const SuffixTree &tree, vector <unsigned int> &suffix);

vector <unsigned int> buildOddSuffixArray(
    const vector <unsigned int> &evenSuffix,
    const vector <int> &input);

vector <unsigned int> buildOddLcp(const RandomLCPGetter &getter,
    const vector <unsigned int> &oddSuffix,
    const vector <int> &input
);

SuffixTree buildSuffixTreeFromSA(vector <unsigned int> &sa,
    vector <unsigned int> &lcp,
    unsigned int length);

SuffixTree buildOddSuffixTree(const SuffixTree &even,
    const vector<int> &input);

void copyNodeExceptParentAndChildren(const SuffixTree::Node &from,
    SuffixTree::Node &to);

void copySubTree(const SuffixTree &from, SuffixTree &to,
    unsigned int fromStart, unsigned int toStart);

unsigned int appendCopyNode(const SuffixTree &from, SuffixTree &to,
    unsigned int toStart, unsigned int u);


template<class Struct, class ActionType>
void doSomething(IndexedPair<Struct> &merging, ActionType action) {
    action(merging[0]);
    action(merging[1]);
}

template<class Struct, class Compare>
Struct &minimal(IndexedPair<Struct> &merging, Compare comp) {
    if (comp(merging[0]) < comp(merging[1])) {
        return merging[0];
    }
    return merging[1];
}

void mergeNodes(unsigned int first, unsigned int second, unsigned int to,
    SuffixTree &result, SuffixTree &tree1, SuffixTree &tree2,
    const vector <int> &input);

bool findSuffixesDfs(unsigned int v, const SuffixTree &tree,
    vector <unsigned int> &output);

vector <unsigned int> findSuffixes(const SuffixTree &tree);

void findLcaSuffix(unsigned int v, const SuffixTree &merged,
    IndexedPair<const SuffixTree &> trees,
    vector <IndexedPair<unsigned int> > &output);

void trueLengthFillDfs(unsigned int v, const SuffixTree &merged,
    const RandomLCPGetter &getter, IndexedPair<unsigned int> suffix, unsigned int length,
    vector <unsigned int> &output);

vector <unsigned int> computeTrueLength(const SuffixTree &merged,
    IndexedPair<const SuffixTree &> trees, unsigned int inputLength);

template<class T>
T *xorPointers(T *a, T *b, T *c) {
    if (a == b) {
        return c;
    }
    return b;
}

void correctMerge(unsigned int v, unsigned int parentsPlace, SuffixTree &merged,
    IndexedPair<MergeTreesStruct> &updatedTrees,
    const vector <unsigned int> &trueLength, const vector <int> &input,
    unsigned int copyTree=2);

void checkTree(const char *message, const SuffixTree &tree,
    const vector <int> &input);

SuffixTree mergeTrees(SuffixTree &tree1, SuffixTree &tree2,
    const vector <int> &input);

int cleanTreeDfs(unsigned int v, unsigned int parent, SuffixTree &tree);

SuffixTree buildTempSuffixTree(const vector <int> &input);

SuffixTree buildSuffixTree(const vector <int> &input);

unsigned long long countSubstrings(SuffixTree &tree, unsigned int v);

void randGen(int length, int alph, int tests);

int main();


#endif
