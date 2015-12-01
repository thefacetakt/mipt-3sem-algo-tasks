// #define _GLIBCXX_DEBUG
// #define _PRINT_DBG

#include <cstdio>
#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>
#include <climits>
#include <cassert>
#include <functional>
#include <string>
#include "LCA.hpp"

using std::vector;
using std::make_pair;
using std::max;
using std::min;
using std::pair;
using std::function;
using std::swap;

#define firstHiddenInfo indexOfParentEdge
#define secondHiddenInfo lastIndex


struct SuffixTree {
    struct Node {
        int parent;
        unsigned int lastIndex;
        unsigned int indexOfParentEdge;
        unsigned int depth;
        int leaf;
        vector <unsigned int> children;
    };
    
    vector <Node> nodes;
    vector <unsigned int> pull;
    unsigned int root;
    
    void checkNode(const Node &node) {
        assert(node.parent == -1 || nodes[node.parent].depth < node.depth);
        assert(node.leaf == -2 || node.parent == -1 || node.indexOfParentEdge < node.lastIndex);
    }
    
    SuffixTree() {
        nodes.resize(1);
        root = 0;
        Node &rootNode = nodes[root];
        rootNode.parent = -1;
        rootNode.lastIndex = 0;
        rootNode.indexOfParentEdge = 0;
        rootNode.leaf = 0;
        rootNode.depth = 0;
    }
    
    void deleteNode(unsigned int index, vector <unsigned int> &newParentsChildren, unsigned int leaf) {
        Node &node = nodes[index];
        assert(node.children.size() == 1);
        for (auto const &u: node.children) {
            nodes[u].parent = node.parent;
            nodes[u].indexOfParentEdge = leaf + nodes[node.parent].depth;
            nodes[u].lastIndex = leaf + nodes[u].depth;
#ifdef _GLIBCXX_DEBUG
            checkNode(nodes[u]);
#endif
            newParentsChildren.push_back(u);
        }
        nodes[index].children.clear();
        pull.push_back(index);
    }
    
    void deleteUselessNode(unsigned int v, unsigned int inParentIndex, int leaf) {
        unsigned int parent = nodes[v].parent;
        nodes[parent].children[inParentIndex] = -1;
        for (auto const &u: nodes[v].children) {
            nodes[parent].children[inParentIndex] = u;
            nodes[u].parent = parent;
            nodes[u].indexOfParentEdge = leaf + nodes[parent].depth;
            nodes[u].lastIndex = leaf + nodes[u].depth;
#ifdef _GLIBCXX_DEBUG
            checkNode(nodes[u]);
#endif
        }
        nodes[v].children.clear();
        pull.push_back(v);
    }
     
    unsigned int newNode() {
        if (pull.size()) {
            unsigned int returnValue = pull.back();
            pull.pop_back();
            return returnValue;
        }
        nodes.push_back(Node());
        return nodes.size() - 1;
    }
    
    
    unsigned int lengthOfEdge(unsigned int index) {
        if (nodes[index].parent == -1) {
            return 0;
        } else {
            return nodes[index].depth - nodes[nodes[index].parent].depth;
        }
    }
    
    unsigned int splitEdge(unsigned int parent, unsigned int indexOfChild, unsigned int length) {
        unsigned int child = nodes[parent].children[indexOfChild];
        unsigned int newNodeIndex = newNode();
        Node &newNode = nodes[newNodeIndex];
        newNode.parent = parent;
        newNode.leaf = -1;
        newNode.indexOfParentEdge = nodes[child].indexOfParentEdge;
        newNode.lastIndex = newNode.indexOfParentEdge + length;
        newNode.depth = nodes[parent].depth + length;
        nodes[child].indexOfParentEdge = newNode.lastIndex;
        nodes[child].parent = newNodeIndex;
        newNode.children.push_back(child);
        nodes[parent].children[indexOfChild] = newNodeIndex;
#ifdef _GLIBCXX_DEBUG
            checkNode(nodes[child]);
            checkNode(nodes[newNodeIndex]);
#endif
        return newNodeIndex;
    }
};

void swapTrees(SuffixTree &first, SuffixTree &second) {
    swap(first.root, second.root);
    first.nodes.swap(second.nodes);
    first.pull.swap(second.pull);
}

SuffixTree buildSuffixTree(const vector<int> &);

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

struct NumberedPair {
    pair<int, int> elements;
    unsigned int number;
    NumberedPair(int first, int second, unsigned int number) : 
        elements(first, second),
        number(number) {
    }
    NumberedPair() {
    }
};

int firstElementOfNumberedPair(const NumberedPair &x) {
   return x.elements.first;
}

int secondElementOfNumberedPair(const NumberedPair &x) {
    return x.elements.second;
}

void depthDfs(unsigned int v, SuffixTree &tree) {
    auto &currentNode = tree.nodes[v];
    if (v != tree.root) {
        currentNode.depth = tree.nodes[currentNode.parent].depth + currentNode.lastIndex - currentNode.indexOfParentEdge;
    }
    for (auto const &u: currentNode.children) {
        depthDfs(u, tree);
    }
}

int decompressDfs(unsigned int v, vector <unsigned int> &parentsNewChildren,
                   SuffixTree &compressed, const vector <int> &input, unsigned int depth) {
    int leaf = -1;
    vector <unsigned int> myNewChildren;
    compressed.nodes[v].depth = depth;
    for (auto const &u: compressed.nodes[v].children) {
        int newLeaf = decompressDfs(u, myNewChildren, compressed, input, min(static_cast<unsigned int>(input.size()), 
                                    depth + (compressed.nodes[u].lastIndex - compressed.nodes[u].indexOfParentEdge) * 2 
                                    - (compressed.nodes[u].lastIndex == input.size() / 2 + 1))); 
        if (leaf == -1) {
            leaf = newLeaf;
        }
    }
    compressed.nodes[v].children = myNewChildren;
    myNewChildren.clear();
    
    compressed.nodes[v].lastIndex = min(compressed.nodes[v].lastIndex * 2, static_cast<unsigned int>(input.size()));
    compressed.nodes[v].indexOfParentEdge *= 2;
    if (compressed.nodes[v].leaf != -1) {
        compressed.nodes[v].leaf = min(compressed.nodes[v].leaf * 2, static_cast<int>(input.size()));
    }
    
    if (compressed.nodes[v].children.size()) {
        vector <unsigned int> similarKids;
        similarKids.push_back(compressed.nodes[v].children[0]);
        for (unsigned int i = 1; i <= compressed.nodes[v].children.size(); ++i) {
            if (i != compressed.nodes[v].children.size()) {
                unsigned int current = compressed.nodes[v].children[i];
                unsigned int previous = compressed.nodes[v].children[i - 1];
                while (i < compressed.nodes[v].children.size() && 
                    input[compressed.nodes[current].indexOfParentEdge] ==
                    input[compressed.nodes[previous].indexOfParentEdge]) {
                    similarKids.push_back(current);
                    ++i;
                    previous = current;
                    if (i < compressed.nodes[v].children.size()) {
                        current = compressed.nodes[v].children[i];
                    }
                }
            }
            if (similarKids.size() != 1) {
                unsigned int newNodeIndex = compressed.newNode();
                auto &newNode = compressed.nodes[newNodeIndex];
                newNode.parent = v;
                newNode.indexOfParentEdge = compressed.nodes[similarKids.back()].indexOfParentEdge;
                newNode.depth = compressed.nodes[v].depth + 1;
                newNode.lastIndex = newNode.indexOfParentEdge + 1;
                newNode.children = similarKids;
                newNode.leaf = -1;
                for (auto &u: similarKids) {
                    auto &currentChild = compressed.nodes[u];
                    currentChild.parent = newNodeIndex;
                    ++currentChild.indexOfParentEdge;
                }
                if (newNode.children.size() && compressed.nodes[newNode.children[0]].lastIndex 
                    == compressed.nodes[newNode.children[0]].indexOfParentEdge) {
                    newNode.leaf = compressed.nodes[newNode.children[0]].leaf;
                    compressed.pull.push_back(newNode.children[0]);
                    newNode.children = vector<unsigned int>(newNode.children.begin() + 1, newNode.children.end());
                }
                myNewChildren.push_back(newNodeIndex);
                compressed.checkNode(compressed.nodes[newNodeIndex]);
                compressed.checkNode(compressed.nodes[v]);
            } else {
                myNewChildren.push_back(similarKids.back());
            }
            if (i != compressed.nodes[v].children.size()) {
                similarKids.clear();
                similarKids.push_back(compressed.nodes[v].children[i]);
            }
        }
    }
    compressed.nodes[v].children = myNewChildren;
    if (myNewChildren.size() == 1 && v != compressed.root && compressed.nodes[v].leaf == -1) {
        compressed.deleteNode(v, parentsNewChildren, leaf);
    } else {
        parentsNewChildren.push_back(v);
    }
    if (compressed.nodes[v].leaf != -1) {
        leaf = compressed.nodes[v].leaf;
    }
    return leaf;
}

void decompress(SuffixTree &compressed, const vector <int> &input) {
    vector <unsigned int> tmp;
    decompressDfs(compressed.root, tmp, compressed, input, 0);
    assert(tmp.size() == 1 && tmp[0] == compressed.root);
    depthDfs(compressed.root, compressed);
}

bool suffixArrayDfs(unsigned int v, unsigned int currentLca,
                    SuffixTree &tree, vector <unsigned int> &array, vector <unsigned int> &lca) {
    auto &currentNode = tree.nodes[v];
    bool returnValue = false;
    if (currentNode.leaf != -1) {
        if (array.size()) {
            lca.push_back(currentLca);
        }
        array.push_back(currentNode.leaf);
        currentLca = currentNode.depth;
        returnValue = true;
    }
    
    for (auto const &u: currentNode.children) {
        if (suffixArrayDfs(u, currentLca, tree, array, lca)) {
            currentLca = currentNode.depth;
            returnValue = true;
        }
    }
    return returnValue;
}

void buildSuffixArray(SuffixTree &tree, vector <unsigned int> &array, vector <unsigned int> &lca) {
    suffixArrayDfs(tree.root, 0, tree, array, lca);
}

void buildEulerDfs(int v, int &count, int depth, const SuffixTree &tree, vector <pair<unsigned int, unsigned int> > &euler, vector <unsigned int> &realNode) {
    auto &currentNode = tree.nodes[v];
    realNode.push_back(v);
    euler.push_back(make_pair(depth, count));
    int currentCount = count;
    for (auto const &u: currentNode.children) {
        buildEulerDfs(u, ++count, depth + 1, tree, euler, realNode);
        euler.push_back(make_pair(depth, currentCount));
    }
}

void buildEuler(const SuffixTree &tree, vector <pair<unsigned int, unsigned int> > &euler, vector <unsigned int> &realNode) {
    int count = 0;
    buildEulerDfs(tree.root, count, 0, tree, euler, realNode);
}

SuffixTree buildSuffixTreeFromSA(vector <unsigned int> &sa, vector <unsigned int> &lcp, unsigned int length) {
    vector <int> tmp;
    SuffixTree result = buildSuffixTree(tmp);
    result.nodes[result.root].leaf = -1;
    int newNodeIndex = result.newNode();
    auto &newNode = result.nodes[newNodeIndex];
    newNode.parent = result.root;
    newNode.indexOfParentEdge = sa[0];
    newNode.lastIndex = length;
    newNode.depth = newNode.lastIndex - newNode.indexOfParentEdge;
    newNode.leaf = sa[0];
    result.nodes[result.root].children.push_back(newNodeIndex);
    unsigned int current = newNodeIndex;
    for (unsigned int i = 1; i < sa.size(); ++i) {
        
        while (result.nodes[current].parent != -1 && result.nodes[result.nodes[current].parent].depth >= lcp[i - 1]) {
            current = result.nodes[current].parent;
        }
        unsigned int parent;
        if (result.nodes[current].parent != -1 && result.nodes[result.nodes[current].parent].depth == lcp[i - 1]) {
            parent = result.nodes[current].parent;
        } else if (result.nodes[current].depth == lcp[i - 1]) {
            parent = current;
        } else {
            parent = result.newNode();
            auto &parentNode = result.nodes[parent];
            parentNode.parent = result.nodes[current].parent;
            parentNode.indexOfParentEdge = result.nodes[current].indexOfParentEdge;
            parentNode.depth = lcp[i - 1];
            parentNode.lastIndex = parentNode.depth - result.nodes[parentNode.parent].depth + parentNode.indexOfParentEdge;
            parentNode.leaf = (length - sa[i] == lcp[i - 1] ? sa[i] : -1);
            result.nodes[result.nodes[current].parent].children.pop_back();
            result.nodes[result.nodes[current].parent].children.push_back(parent);
            parentNode.children.push_back(current);
            result.nodes[current].parent = parent;
            result.nodes[current].indexOfParentEdge += parentNode.depth - result.nodes[parentNode.parent].depth;
        }
        if (lcp[i - 1] != length - sa[i]) {
            newNodeIndex = result.newNode();
            result.nodes[parent].children.push_back(newNodeIndex);
            auto &newNode1 = result.nodes[newNodeIndex];
            newNode1.parent = parent;
            newNode1.leaf = sa[i];
            newNode1.lastIndex = length;
            newNode1.depth = length - sa[i];
            newNode1.indexOfParentEdge = sa[i] + lcp[i - 1];
            parent = newNodeIndex;
        }
        current = parent;
    }
    return result;
}

SuffixTree buildOddSuffixTree(SuffixTree &even, const vector <int> &input) {
    vector <unsigned int> evenSuffix;
    vector <unsigned int> evenLca;
    
    buildSuffixArray(even, evenSuffix, evenLca);
    vector <pair<unsigned int, unsigned int> > evenEuler;
    vector <unsigned int> evenRealNode;
    buildEuler(even, evenEuler, evenRealNode);
    
    LCA evenRandomLCA(evenEuler);
    
    vector <unsigned int> evenIrrealNode(input.size());
    for (unsigned int i = 0; i < evenRealNode.size(); ++i) {
        int leaf = even.nodes[evenRealNode[i]].leaf;
        if (leaf != -1)
            evenIrrealNode[leaf] = i;
    }

    vector <unsigned int> antiSuffixArray(input.size());
    for (unsigned int i = 0; i < evenSuffix.size(); ++i) {
        antiSuffixArray[evenSuffix[i]] = i;
    }
    vector <unsigned int> oddSuffix;
    for (unsigned int i = 1; i < input.size(); i += 2) {
        oddSuffix.push_back(i);
    }
    
    oddSuffix = radixSort(radixSort(oddSuffix, [&input, &antiSuffixArray] (unsigned int i) -> int {
        return (i + 1 == input.size() ? -1 : antiSuffixArray[i + 1]);
    }), [&input] (unsigned int i) -> unsigned int {
        return input[i];
    });
    
    vector <unsigned int> oddLca(oddSuffix.size() - 1);
    for (unsigned int i = 0; i + 1 < oddSuffix.size(); ++i) {
        if (input[oddSuffix[i]] == input[oddSuffix[i + 1]]) {
            if (oddSuffix[i + 1] + 1 == input.size() || oddSuffix[i] + 1 == input.size()) {
                oddLca[i] = 1;
            } else {
                oddLca[i] = even.nodes[evenRealNode[evenRandomLCA.lca(evenIrrealNode[oddSuffix[i] + 1], evenIrrealNode[oddSuffix[i + 1] + 1])]].depth + 1;
            }
        } else {
            oddLca[i] = 0;
        }
    }
#ifdef _PRINT_DBG
    for (unsigned int i = 0; i < oddSuffix.size(); ++i) {
        if (i + 1 != oddSuffix.size()) {
            assert(vector<int>(input.begin() + oddSuffix[i], input.begin() + oddSuffix[i] + oddLca[i])
                    == vector<int>(input.begin() + oddSuffix[i + 1], input.begin() + oddSuffix[i + 1] + oddLca[i])
                    &&
                    (oddSuffix[i] + oddLca[i] == input.size() || oddSuffix[i + 1] + oddLca[i] == input.size() || 
                vector<int>(input.begin() + oddSuffix[i], input.begin() + oddSuffix[i] + oddLca[i] + 1)
                != vector<int>(input.begin() + oddSuffix[i + 1], input.begin() + oddSuffix[i + 1] + oddLca[i] + 1)));
        }
        for (unsigned int j = oddSuffix[i]; j < input.size(); ++j) {
            printf("%d ", input[j]);
        }
        printf("| %d\n", (i == oddLca.size() ? -1 : oddLca[i]));
    }
#endif
    return buildSuffixTreeFromSA(oddSuffix, oddLca, input.size());
}

inline void copyNodeExceptParentAndChildren(const SuffixTree::Node &from,SuffixTree::Node &to) {
    to.lastIndex = from.lastIndex;
    to.depth = from.depth;
    to.indexOfParentEdge = from.indexOfParentEdge;
    to.leaf = from.leaf;
}

void copySubTree(const SuffixTree &from, SuffixTree &to, unsigned int fromStart, unsigned int toStart);

unsigned int appendCopyNode(const SuffixTree &from, SuffixTree &to, unsigned int toStart, unsigned int u) {
    unsigned int newStart = to.newNode();
    to.nodes[newStart].parent = toStart;
    to.nodes[toStart].children.push_back(newStart);
    copyNodeExceptParentAndChildren(from.nodes[u], to.nodes[newStart]);
    copySubTree(from, to, u, newStart);
    return newStart;
}

void copySubTree(const SuffixTree &from, SuffixTree &to, unsigned int fromStart, unsigned int toStart) {
    for (auto const &u: from.nodes[fromStart].children) {
        appendCopyNode(from, to, toStart, u);
    }
}

void mergeNodes(unsigned int first, unsigned int second, unsigned int to,
    SuffixTree &result, SuffixTree &tree1, SuffixTree &tree2, const vector <int> &input) {
    unsigned int firstIndex = 0;
    unsigned int secondIndex = 0;
    for (; firstIndex < tree1.nodes[first].children.size() && secondIndex < tree2.nodes[second].children.size();) {
        unsigned int firstChild = tree1.nodes[first].children[firstIndex];
        unsigned int secondChild = tree2.nodes[second].children[secondIndex];
        int firstChar = input[tree1.nodes[firstChild].indexOfParentEdge];
        int secondChar = input[tree2.nodes[secondChild].indexOfParentEdge];
        if (firstChar < secondChar) {
            appendCopyNode(tree1, result, to, firstChild);
            ++firstIndex;
            continue;
        } else if (firstChar > secondChar) {
            appendCopyNode(tree2, result, to, secondChild);
            ++secondIndex;
            continue;
        } else {
            int firstLength = tree1.lengthOfEdge(firstChild);
            int secondLength = tree2.lengthOfEdge(secondChild);
            
            if (firstLength < secondLength) {
                secondChild = tree2.splitEdge(second, secondIndex, firstLength);
                
            } else if (firstLength > secondLength) {
                firstChild = tree1.splitEdge(first, firstIndex, secondLength);
            }
            unsigned int newNodeIndex = result.newNode();
            
            auto &newNode = result.nodes[newNodeIndex];
            newNode.parent = to;
            result.nodes[to].children.push_back(newNodeIndex);
            if (min(firstLength, secondLength) == 1) {
                if (firstLength < secondLength || (firstLength == secondLength && tree1.nodes[firstChild].leaf != -1)) {
                    copyNodeExceptParentAndChildren(tree1.nodes[firstChild], newNode);
                } else {
                    copyNodeExceptParentAndChildren(tree2.nodes[secondChild], newNode);
                } 
            } else {
                newNode.leaf = -2;
                newNode.depth = min(firstLength, secondLength) + result.nodes[to].depth;
                newNode.firstHiddenInfo = firstChild;
                newNode.secondHiddenInfo = secondChild;
            }
            mergeNodes(firstChild, secondChild, newNodeIndex, result, tree1, tree2, input);
            ++firstIndex;
            ++secondIndex;
            continue;
        }
    }
    for (; firstIndex < tree1.nodes[first].children.size(); ++firstIndex) {
        appendCopyNode(tree1, result, to, tree1.nodes[first].children[firstIndex]);
    }
    for (; secondIndex < tree2.nodes[second].children.size(); ++secondIndex) {
        appendCopyNode(tree2, result, to, tree2.nodes[second].children[secondIndex]);
    }
}


bool findSuffixesDfs(unsigned int v, SuffixTree &tree, vector <unsigned int> &output) {
    while (output.size() <= v)
        output.push_back(-1);
    output[v] = tree.nodes[v].leaf;
    for (auto const &u: tree.nodes[v].children) {
        if (findSuffixesDfs(u, tree, output)) {
            if (output[v] == -1) {
                output[v] = output[u];
            }
        }
    }
    return (output[v] != -1);
}

vector <unsigned int> findSuffixes(SuffixTree &tree) {
    vector <unsigned int> result;
    findSuffixesDfs(tree.root, tree, result);
    return result;
}

void findLcaSuffix(unsigned int v, const SuffixTree &merged, const SuffixTree &first, const SuffixTree &second, vector <pair<int, int> > &output) {
    while (output.size() <= v) {
        output.push_back(make_pair(-1, -1));
    }
    if (merged.nodes[v].leaf >= 0) {
        if (merged.nodes[v].leaf % 2 == 0) {
            output[v].first = merged.nodes[v].leaf;
        } else {
            output[v].second = merged.nodes[v].leaf;
        }
    } else if (merged.nodes[v].leaf == -2) {
        output[v].first = first.nodes[merged.nodes[v].firstHiddenInfo].leaf;
        output[v].second = second.nodes[merged.nodes[v].secondHiddenInfo].leaf;
    }
    int reserve = -1;
    for (auto const &u: merged.nodes[v].children) {
        findLcaSuffix(u, merged, first, second, output);
        if (output[u].first != -1) {
            if (output[v].first == -1) {
                output[v].first = output[u].first;
                if (output[u].second != -1) {
                    reserve = output[u].second;
                }
            } else if (reserve != -1) {
                output[v].first = output[u].first;
                output[v].second = reserve;
            }
        }
        if (output[u].second != -1 && output[u].second != reserve) {
            if (output[v].second == -1) {
                output[v].second = output[u].second;
            }
        }
    }
    assert(merged.nodes[v].leaf != -2 || (output[v].first != -1 && output[v].second != -1));
}

void trueLengthFillDfs(unsigned int v, const SuffixTree &merged, 
                       const function<unsigned int(unsigned int, unsigned int)> &suffixLcaGetter,
                       unsigned int firstSuffix, unsigned int secondSuffix, unsigned int length,
                       vector <unsigned int> &output) {
    while (output.size() <= v) {
        output.push_back(-1);
    }
    if (output[v] != -1) {
        return;
    }
    if (v == merged.root) {
        output[v] = 0;
        return;
    }
    if (firstSuffix + 1 == length || secondSuffix + 1 == length) {
        output[v] = 1;
        return;
    }
    
    unsigned int u = suffixLcaGetter(firstSuffix + 1, secondSuffix + 1);
    trueLengthFillDfs(suffixLcaGetter(firstSuffix + 1, secondSuffix + 1), 
                        merged, suffixLcaGetter, firstSuffix + 1, 
                        secondSuffix + 1, length, output);
    output[v] = output[u] + 1;
    assert(output[v] <= merged.nodes[v].depth);
}

vector <unsigned int> computeTrueLength(const SuffixTree &merged, const SuffixTree &first, const SuffixTree &second,
                       unsigned int inputLength) {
    vector <pair<unsigned int, unsigned int> > euler;
    vector <unsigned int> realNode;
    buildEuler(merged, euler, realNode);
    
    vector <unsigned int> irrealNode(inputLength);
    for (unsigned int i = 0; i < realNode.size(); ++i) {
        int leaf = merged.nodes[realNode[i]].leaf;
        if (leaf == -2) {
            int leaf1 = first.nodes[merged.nodes[realNode[i]].firstHiddenInfo].leaf;
            if (leaf1 != -1) {
                irrealNode[leaf1] = i;
            }
            leaf = second.nodes[merged.nodes[realNode[i]].secondHiddenInfo].leaf;
        } 
        if (leaf != -1) {
            irrealNode[leaf] = i;
        }
    }
    
    LCA lcaGetter(euler);
    
    vector <unsigned int> output;
    
    function<unsigned int(unsigned int, unsigned int)> suffixLcaGetter = [&lcaGetter, &realNode, &irrealNode] (unsigned int i, unsigned int j) -> unsigned int {
        return realNode[lcaGetter.lca(irrealNode[i], irrealNode[j])];
    };
    
    vector <pair<int, int> > superSuffix;
    findLcaSuffix(merged.root, merged, first, second, superSuffix);
    
    vector <unsigned int> toVisit;
    toVisit.push_back(merged.root);
    while (toVisit.size()) {
        int v = toVisit.back();
        toVisit.pop_back();
        if (merged.nodes[v].leaf == -2 && (output.size() <= v || output[v] == -1)) {
            trueLengthFillDfs(v, merged, suffixLcaGetter, 
                              superSuffix[v].first,
                              superSuffix[v].second,
                              inputLength, output);
        }
        for (auto const &u: merged.nodes[v].children) {
            toVisit.push_back(u);
        }
    }
    return output;
}

void correctMerge(unsigned int v, unsigned int parentsPlace, SuffixTree &merged, const SuffixTree &first, 
                  const SuffixTree &second, const vector <unsigned int> &trueLength, const vector <int> &input, 
                  const vector <unsigned int> &firstSuffix, const vector <unsigned int> &secondSuffix,
                  unsigned int copyTree=2) {
    if (merged.nodes[v].leaf == -2) {
        unsigned int firstInfo = merged.nodes[v].firstHiddenInfo;
        unsigned int secondInfo = merged.nodes[v].secondHiddenInfo;
        SuffixTree::Node const *firstNode = &first.nodes[firstInfo];
        SuffixTree::Node const *secondNode = &second.nodes[secondInfo];
        if (copyTree == 0 || (merged.nodes[v].depth == trueLength[v] && firstNode->leaf != -1)) {
            copyNodeExceptParentAndChildren(*firstNode, merged.nodes[v]);
        } else if (copyTree == 1 || merged.nodes[v].depth == trueLength[v]) {
            copyNodeExceptParentAndChildren(*secondNode, merged.nodes[v]);
        } else {
            unsigned int commonLength = trueLength[v] - merged.nodes[merged.nodes[v].parent].depth;
            if (input[firstSuffix[firstInfo] + (firstNode->parent == -1 ? 0 : first.nodes[firstNode->parent].depth) + commonLength] >
                input[secondSuffix[secondInfo] + (secondNode->parent == -1 ? 0 : second.nodes[secondNode->parent].depth) + commonLength]) {
                
                swap(firstInfo, secondInfo);
                swap(firstNode, secondNode);
                copyTree = 1;
            } else {
                copyTree = 0;
            }
            copyNodeExceptParentAndChildren(*firstNode, merged.nodes[v]);
            unsigned int newNodeIndex = merged.splitEdge(merged.nodes[v].parent, parentsPlace, commonLength);
            unsigned int newCopy = merged.newNode();
            merged.nodes[newCopy].parent = newNodeIndex;
            merged.nodes[newNodeIndex].children.push_back(newCopy);
            copyNodeExceptParentAndChildren(*secondNode, merged.nodes[newCopy]);
            merged.nodes[newCopy].indexOfParentEdge += commonLength;
            copySubTree((copyTree == 0 ? second : first), merged, secondInfo, newCopy);
        }
    }
    if (copyTree != 2 && merged.nodes[v].leaf != -1 && merged.nodes[v].leaf % 2 != copyTree) {
        merged.nodes[v].leaf = -1;
    }
    for (unsigned int i = 0; i < merged.nodes[v].children.size(); ++i) {
        unsigned int u = merged.nodes[v].children[i];
        correctMerge(u, i, merged, first, second, trueLength, input, firstSuffix, secondSuffix, copyTree);
    }
}

void dfs__(SuffixTree &tree, int v) {
#ifdef _PRINT_DBG
    printf("%d -> %d [label=\"%d:%d, %d %d\"]\n", tree.nodes[v].parent, v, tree.nodes[v].indexOfParentEdge, tree.nodes[v].lastIndex, tree.nodes[v].leaf, tree.nodes[v].depth);
#endif
    for (auto const &u: tree.nodes[v].children) {
        dfs__(tree, u);
    }
}


void dfs_(SuffixTree &tree, int v, int count, const vector<int> &input) {
#ifdef _PRINT_DBG
    printf("%d -> %d [label=\"", tree.nodes[v].parent, v);
    for (int i = tree.nodes[v].indexOfParentEdge; i != tree.nodes[v].lastIndex; ++i) {
        printf("%d", input[i]);
    }
    printf(", %d, %d\"]\n", tree.nodes[v].leaf, tree.nodes[v].depth);
#endif
    assert(tree.nodes[v].indexOfParentEdge < tree.nodes[v].lastIndex || tree.nodes[v].parent == -1);
    assert(tree.nodes[v].parent == -1 || tree.nodes[v].depth > tree.nodes[tree.nodes[v].parent].depth);
    for (auto const &u: tree.nodes[v].children) {
        dfs_(tree, u, count + 1, input);
    }
}


SuffixTree mergeTrees(SuffixTree &tree1, SuffixTree &tree2, const vector <int> &input) {
    vector <int> tmp;
    SuffixTree merged = buildSuffixTree(tmp);
    merged.nodes[merged.root].leaf = -1;
    mergeNodes(tree1.root, tree2.root, merged.root, merged, tree1, tree2, input);

#ifdef _PRINT_DBG
    dfs__(merged, merged.root);
    printf("NEW TREE1\n");
    dfs_(tree1, tree1.root, 0, input);
    printf("NEW TREE2\n");
    dfs_(tree2, tree2.root, 0, input);
#endif

    vector <unsigned int> firstSuffixes = findSuffixes(tree1);
    vector <unsigned int> secondSuffixes = findSuffixes(tree2);
    
    vector <unsigned int> trueLength = computeTrueLength(merged, tree1, tree2, input.size());
    correctMerge(merged.root, -1, merged, tree1, tree2, trueLength, input, firstSuffixes, secondSuffixes);
    
    return merged;
}

int cleanTreeDfs(unsigned int v, unsigned int parent, SuffixTree &tree) {
    int leaf = tree.nodes[v].leaf;
    for (unsigned int i = 0; i < tree.nodes[v].children.size(); ++i) {
        unsigned int u = tree.nodes[v].children[i];
        int newLeaf = cleanTreeDfs(u, i, tree);
        if (leaf == -1) {
            leaf = newLeaf;
        }
    }
    int freeCell = 0;
    for (unsigned int i = 0; i < tree.nodes[v].children.size(); ++i) {
        if (tree.nodes[v].children[i] != -1) {
            tree.nodes[v].children[freeCell] = tree.nodes[v].children[i];
            ++freeCell;
        }
    }
    tree.nodes[v].children.resize(freeCell);
    if (tree.nodes[v].children.size() <= 1 && tree.nodes[v].leaf == -1 && v != tree.root) {
        tree.deleteUselessNode(v, parent, leaf);
    }
    return leaf;
}



SuffixTree buildSuffixTree(const vector <int> &input) {
    if (input.size() == 0) {
        return SuffixTree();
    }
    if (input.size() == 1) {
        SuffixTree result;
        unsigned int newNodeIndex = result.newNode();
        SuffixTree::Node &newNode = result.nodes[newNodeIndex];
        newNode.parent = result.root;
        newNode.lastIndex = 1;
        newNode.indexOfParentEdge = 0;
        newNode.depth = 1;
        newNode.leaf = 0;
        result.nodes[result.root].children.push_back(newNodeIndex);
        result.nodes[result.root].leaf = 1;
        return result;
    }
    vector <NumberedPair> consequentPairs;
    for (unsigned int i = 0; i < input.size(); i += 2) {
        consequentPairs.push_back(NumberedPair(
                input[i], 
                (i + 1 == input.size() ? -1: input[i + 1]),
                i / 2
            )
        );        
    }
    consequentPairs = radixSort(radixSort(consequentPairs, secondElementOfNumberedPair), firstElementOfNumberedPair);
    vector <int> pairedString((input.size() + 1) / 2);
    unsigned int currentNumber = 0;
    for (unsigned int i = 0; i < consequentPairs.size(); ++i) {
        if (i && consequentPairs[i].elements != consequentPairs[i - 1].elements) {
            ++currentNumber;
        }
        pairedString[consequentPairs[i].number] = currentNumber;
    }
#ifdef _PRINT_DBG
    for (unsigned int i = 0; i < pairedString.size(); ++i) {
        printf("%d ", pairedString[i]);
    }
    printf("\n");
#endif
    SuffixTree compressed = buildSuffixTree(pairedString);
    compressed.nodes[compressed.root].leaf = -1;
    decompress(compressed, input);
    SuffixTree &even = compressed;
#ifdef _PRINT_DBG
    printf("EVEN:\n");
    dfs_(even, even.root, 0, input);
#endif
    SuffixTree odd = buildOddSuffixTree(even, input);
#ifdef _PRINT_DBG
    printf("ODD:\n");
    dfs_(odd, odd.root, 0, input);
#endif
    SuffixTree almostResult = mergeTrees(even, odd, input);
#ifdef _PRINT_DBG
    printf("ALMOST:\n");
    dfs_(almostResult, almostResult.root, 0, input);
#endif
    cleanTreeDfs(almostResult.root, -1, almostResult);
    almostResult.nodes[almostResult.root].leaf = input.size();
#ifdef _PRINT_DBG
    printf("RESL\n");
    dfs_(almostResult, almostResult.root, 0, input);
#endif
    return almostResult;
}

long long dfs(SuffixTree &t, int v) {
    if (v != t.root)
        assert(t.nodes[v].depth > t.nodes[t.nodes[v].parent].depth);
    long long sum = (v == t.root ? 0 : t.nodes[v].depth - t.nodes[t.nodes[v].parent].depth);
    for (auto const &u: t.nodes[v].children) {
        sum += dfs(t, u);
    }
    return sum;
}

vector <int> res;
void gen(int i, int n) {
    if (i == n) {
        for (int j = 0; j < n; ++j) {
            printf("%d ", res[j]);
        }
        printf("\n");
        buildSuffixTree(res);
        return;
    }
    for (int j = 0; j <= 1; ++j) {
        res.push_back(j);
        gen(i + 1, n);
        res.pop_back();
    }
}

int main() {
    std::string s;
    std::cin >> s;
    
    int n = s.size();
    vector <int> x(n);
    for (int i = 0; i < n; ++i)
        x[i] = s[i] - 'a';
    SuffixTree y = buildSuffixTree(x);
    std::cout << dfs(y, y.root) << std::endl;
    return 0;
    
}