#ifndef _USEFUL_STRUCTURES
#define _USEFUL_STRUCTURES

class RandomLCPGetter;
struct NumberedPair;

template<class T>
class IndexedPair;

struct MergeNodesStruct;
struct MergeTreesStruct;

#include "SuffixTree.hpp"
#include "farach.hpp"
#include <vector>

using std::vector;



enum ELeafGetterMode {
    NUMBER,
    LEAF,
};

class RandomLCPGetter {
private:
    const SuffixTree &tree;
    vector <EulerPair> euler;
    vector <unsigned int> realNode;
    vector <unsigned int> irrealNode;
    LCA lcaGetter;
public:
    template<class LeafGetter>
    RandomLCPGetter(const SuffixTree &tree, unsigned int inputLength,
        LeafGetter leafGetter
    ) :
        tree(tree)
    {
        buildEuler(tree, euler, realNode);
        irrealNode.resize(inputLength + 1);
        for (unsigned int i = 0; i < realNode.size(); ++i) {
            for (int j = 0; j < leafGetter(NUMBER); ++j) {
                int leaf = leafGetter(LEAF, realNode[i], j);
                    // tree[realNode[i]].leaf;
                if (leaf != -1) {
                    irrealNode[leaf] = i;
                }
            }
        }
        lcaGetter = LCA(euler);
    }

    unsigned int lca(unsigned int i, unsigned int j) const;

    unsigned int lcp(unsigned int i, unsigned int j) const;
};

struct NumberedPair {
    pair<int, int> elements;
    unsigned int number;
    NumberedPair(int first, int second, unsigned int number);

    NumberedPair();
};

template<class T>
class IndexedPair {
    T first;
    T second;
public:
    IndexedPair(T first, T second) : first(first), second(second) {
    }

    T &operator[](size_t i) {
        #ifdef _GLIBCXX_DEBUG
            assert(0 <= i && i <= 1);
        #endif
        return (i == 0 ? first : second);
    }

    const T &operator[](size_t i) const {
        #ifdef _GLIBCXX_DEBUG
            assert(0 <= i && i <= 1);
        #endif
        return (i == 0 ? first : second);
    }
};

struct MergeNodesStruct {
    SuffixTree &tree;
    const vector <int> &input;
    unsigned int v;
    unsigned int childIndex;
    unsigned int child;
    unsigned int length;
    int symbol;

    MergeNodesStruct(SuffixTree &tree,
        unsigned int v, const vector <int> &input);

    bool end() const;

    void evaluate();

    void split(unsigned int splitLength);
};

struct MergeTreesStruct {
    const SuffixTree &tree;
    vector <unsigned int> suffix;
    unsigned int number;
    unsigned int info;

    MergeTreesStruct(const SuffixTree &tree, unsigned int number);

    void evaluate(const SuffixTree &merged, unsigned int v);
};


#endif
