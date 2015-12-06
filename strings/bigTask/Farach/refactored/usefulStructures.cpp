#include "usefulStructures.hpp"
#include "SuffixTree.hpp"
#include "farach.hpp"


unsigned int RandomLCPGetter::lca(unsigned int i, unsigned int j) const {
    return realNode[lcaGetter.lca(irrealNode[i], irrealNode[j])];
}

unsigned int RandomLCPGetter::lcp(unsigned int i, unsigned int j) const {
    return tree[lca(i, j)].depth;
}


NumberedPair::NumberedPair(int first, int second, unsigned int number) :
    elements(first, second),
    number(number) {
}

NumberedPair::NumberedPair() {
}

MergeNodesStruct::MergeNodesStruct(SuffixTree &tree,
    unsigned int v, const vector <int> &input
) : tree(tree), v(v), input(input) {
    childIndex = 0;
}

bool MergeNodesStruct::end() const {
    return childIndex == tree[v].size();
}

void MergeNodesStruct::evaluate() {
    child = tree[v][childIndex];
    symbol = input[tree[child].indexOfParentEdge];
}

void MergeNodesStruct::split(unsigned int splitLength) {
    child = tree.splitEdge(v, childIndex, splitLength);
}


MergeTreesStruct::MergeTreesStruct(const SuffixTree &tree, unsigned int number)
    : tree(tree), number(number)
{
    suffix = findSuffixes(tree);
}

void MergeTreesStruct::evaluate(const SuffixTree &merged, unsigned int v) {
    info = merged[v].getHiddenInfo(number);
}
