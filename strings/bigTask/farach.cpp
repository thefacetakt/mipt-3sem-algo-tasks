#define _GLIBCXX_DEBUG

#include <cstdio>
#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>
#include <climits>
#include <cassert>
#include "LCA.hpp"

using std::vector;
using std::make_pair;
using std::max;
using std::min;
using std::pair;

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
    
    void deleteNode(unsigned int index, vector <unsigned int> &newParentsChildren) {
        Node &node = nodes[index];
        for (auto const &u: node.children) {
            nodes[u].parent = node.parent;
            nodes[u].indexOfParentEdge = node.indexOfParentEdge;
            newParentsChildren.push_back(u);
        }
        pull.push_back(index);
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
        return newNodeIndex;
    }
};

SuffixTree buildSuffixTree(const vector<int> &);

template<class T, class Compare>
vector <T> radixSort(const vector <T> &input, Compare comp) {
    int maxElement = 0;
    for (auto const &element: input) {
        maxElement = max(maxElement, static_cast<int>(comp(element)));
    }
    vector <unsigned int> count(maxElement + 1);
    for (auto const &element: input) {
        ++count[comp(element)];
    }
    for (unsigned int i = 1; i <= maxElement; ++i) {
        count[i] += count[i - 1];
    }
    vector <T> result(input.size());
    for (unsigned int i = input.size() - 1; i != UINT_MAX; --i) {
        const T &element = input[i];
        result[--count[comp(element)]] = element;
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

void decompressDfs(unsigned int v, vector <unsigned int> &parentsNewChildren,
                   SuffixTree &compressed, const vector <int> &input) {
    vector <unsigned int> myNewChildren;
    
    for (auto const &u: compressed.nodes[v].children) {
        decompressDfs(u, myNewChildren, compressed, input);
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
                    newNode.children = vector<unsigned int>(newNode.children.begin() + 1, newNode.children.end());
                }
                myNewChildren.push_back(newNodeIndex);
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
        compressed.deleteNode(v, parentsNewChildren);
    } else {
        parentsNewChildren.push_back(v);
    }
}

void decompress(SuffixTree &compressed, const vector <int> &input) {
    vector <unsigned int> tmp;
    decompressDfs(compressed.root, tmp, compressed, input);
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

void buildEulerDfs(int v, int &count, int depth, SuffixTree &tree, vector <pair<unsigned int, unsigned int> > &euler, vector <unsigned int> &realNode) {
    auto &currentNode = tree.nodes[v];
    realNode.push_back(v);
    euler.push_back(make_pair(depth, count));
    int currentCount = count;
    for (auto const &u: currentNode.children) {
        buildEulerDfs(u, ++count, depth + 1, tree, euler, realNode);
        euler.push_back(make_pair(depth, currentCount));
    }
}

void buildEuler(SuffixTree &tree, vector <pair<unsigned int, unsigned int> > &euler, vector <unsigned int> &realNode) {
    int count = 0;
    buildEulerDfs(tree.root, count, 0, tree, euler, realNode);
}

SuffixTree buildSuffixTreeFromSA(vector <unsigned int> &sa, vector <unsigned int> &lcp, unsigned int length) {
    vector <int> tmp;
    SuffixTree result = buildSuffixTree(tmp);
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
            parentNode.leaf = (length - sa[i] == lcp[i - 1] ? i : -1);
            result.nodes[result.nodes[current].parent].children.pop_back();
            result.nodes[result.nodes[current].parent].children.push_back(parent);
            parentNode.children.push_back(current);
            result.nodes[current].parent = parent;
            result.nodes[current].indexOfParentEdge = length - (result.nodes[current].depth - parentNode.depth);
        }
        if (lcp[i - 1] != length - sa[i]) {
            newNodeIndex = result.newNode();
            result.nodes[parent].children.push_back(newNodeIndex);
            auto &newNode1 = result.nodes[newNodeIndex];
            newNode1.parent = parent;
            newNode1.leaf = i;
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
    for (int i = 0; i < oddSuffix.size(); ++i) {
        for (int j = oddSuffix[i]; j < input.size(); ++j) {
            printf("%d ", input[j]);
        }
        printf("| %d\n", (i == oddLca.size() ? -1 : oddLca[i]));
    }
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
            newNode.leaf = -2;
            newNode.depth = min(firstLength, secondLength);
            newNode.indexOfParentEdge = firstChild;
            newNode.lastIndex = secondChild;
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

SuffixTree mergeTrees(SuffixTree &tree1, SuffixTree &tree2, const vector <int> &input) {
    vector <int> tmp;
    SuffixTree result = buildSuffixTree(tmp);
    mergeNodes(tree1.root, tree2.root, result.root, result, tree1, tree2, input);
    return result;
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
    for (unsigned int i = 0; i < pairedString.size(); ++i) {
        printf("%d ", pairedString[i]);
    }
    printf("\n");
    SuffixTree compressed = buildSuffixTree(pairedString);
    compressed.nodes[compressed.root].leaf = -1;
    decompress(compressed, input);
    
    SuffixTree &even = compressed;
    SuffixTree odd = buildOddSuffixTree(even, input);
     
//     
//     vector <unsigned int> evenSuffix;
//     vector <unsigned int> evenLca;
//     
//     buildSuffixArray(even, evenSuffix, evenLca);
//     vector <pair<unsigned int, unsigned int> > evenEuler;
//     vector <unsigned int> evenRealNode;
//     buildEuler(tree, evenEuler, evenRealNode);
//     LCA evenRandomLCA(evenEuler);
//     vector <unsigned int> evenIrrealNode(input.size());
//     for (unsigned int i = 0; i < evenReadNode.size(); ++i) {
//         int leaf = even.nodes[evenReadNode[i]];
//         if (leaf != -1)
//             evenIrrealNode[leaf] = i;
//     }
//     for (int i = 0; i < evenIrrealNode.size(); ++i) {
//         printf("%d ", evenIrrealNode[i]);
//     }
//     printf("\n");
    return SuffixTree();
    

    
    return SuffixTree();
}


void dfs_(SuffixTree &tree, int v, int count) {
    for (int i = 0; i < count; ++i) {
        printf("-");
    }
    printf(" %d: %d %d %d %d %d\n", v, tree.nodes[v].parent, tree.nodes[v].indexOfParentEdge, tree.nodes[v].lastIndex, tree.nodes[v].leaf, tree.nodes[v].depth);
    for (auto const &u: tree.nodes[v].children) {
        dfs_(tree, u, count + 1);
    }
}

void sample_() {
    SuffixTree sample;
    sample.nodes[sample.root].leaf = -1;
    int newNode = sample.newNode();
    sample.nodes[newNode].parent = sample.root;
    sample.nodes[newNode].indexOfParentEdge = 1;
    sample.nodes[newNode].lastIndex = 6;
    sample.nodes[newNode].leaf = 1;
    sample.nodes[sample.root].children.push_back(newNode);
    
    
    newNode = sample.newNode();
    sample.nodes[newNode].parent = sample.root;
    sample.nodes[newNode].indexOfParentEdge = 0;
    sample.nodes[newNode].lastIndex = 1;
    sample.nodes[newNode].leaf = -1;
    sample.nodes[sample.root].children.push_back(newNode);
    int root = newNode;
    
    newNode = sample.newNode();
    sample.nodes[newNode].parent = root;
    sample.nodes[newNode].indexOfParentEdge = 1;
    sample.nodes[newNode].lastIndex = 6;
    sample.nodes[newNode].leaf = 0;
    sample.nodes[root].children.push_back(newNode);
    
    newNode = sample.newNode();
    sample.nodes[newNode].parent = root;
    sample.nodes[newNode].indexOfParentEdge = 3;
    sample.nodes[newNode].lastIndex = 6;
    sample.nodes[newNode].leaf = 2;
    sample.nodes[root].children.push_back(newNode);
    
    
    newNode = sample.newNode();
    sample.nodes[newNode].parent = sample.root;
    sample.nodes[newNode].indexOfParentEdge = 5;
    sample.nodes[newNode].lastIndex = 6;
    sample.nodes[newNode].leaf = 5;
    sample.nodes[sample.root].children.push_back(newNode);
    
    newNode = sample.newNode();
    sample.nodes[newNode].parent = sample.root;
    sample.nodes[newNode].indexOfParentEdge = 3;
    sample.nodes[newNode].lastIndex = 6;
    sample.nodes[newNode].leaf = 3;
    sample.nodes[sample.root].children.push_back(newNode);
    
    
    newNode = sample.newNode();
    sample.nodes[newNode].parent = sample.root;
    sample.nodes[newNode].indexOfParentEdge = 4;
    sample.nodes[newNode].lastIndex = 6;
    sample.nodes[newNode].leaf = 4;
    sample.nodes[sample.root].children.push_back(newNode);
    
    vector <int> input;
    input.push_back(1);
    input.push_back(2);
    input.push_back(1);
    input.push_back(1);
    input.push_back(1);
    input.push_back(2);
    input.push_back(2);
    input.push_back(1);
    input.push_back(2);
    input.push_back(2);
    input.push_back(2);
    
    depthDfs(sample.root, sample);
    
    decompress(sample, input);
    
    dfs_(sample, sample.root, 0);
    
    
    SuffixTree &even = sample;
    SuffixTree odd = buildOddSuffixTree(even, input);
    dfs_(odd, odd.root, 0);
    
    SuffixTree merged = mergeTrees(even, odd, input);
    dfs_(even, even.root, 0);
    dfs_(odd, odd.root, 0);
    dfs_(merged, merged.root, 0);
    
    
//     vector <unsigned int> evenSuffixArray;
//     vector <unsigned int> evenLcaArray;
//     buildSuffixArray(sample, evenSuffixArray, evenLcaArray);
//     
//     
//     for (int i = 0; i < evenSuffixArray.size(); ++i) {
//         for (int j = evenSuffixArray[i]; j < input.size(); ++j) {
//             printf("%d ", input[j]);
//         }
//         printf(": %d\n", (i < evenLcaArray.size() ? evenLcaArray[i] : -1));
//     }
}


int main() {
    sample_();
    return 0;
    int n;
    scanf("%d", &n);
    vector <int> x(n);
    for (int i = 0; i < n; ++i)
        scanf("%d", &x[i]);
    buildSuffixTree(x);
    
//     vector <pair<int, int> > x(n);
//     for (int i = 0; i < n; ++i) {
//         scanf("%d %d", &x[i].first, &x[i].second);
//     }
//     x = radixSort(x, [] (const pair<int, int> & y) -> int {
//         return y.second;
//     });
//     printf("---------------\n");
//     for (int i = 0; i < n; ++i) {
//         printf("%d %d\n", x[i].first, x[i].second);
//     }
//     
//     x = radixSort(x, [] (const pair<int, int> & y) -> int {
//         return y.first;
//     });
//     printf("---------------\n");
//     for (int i = 0; i < n; ++i) {
//         printf("%d %d\n", x[i].first, x[i].second);
//     }
    
}