#include <cstdio>
#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>
#include <climits>
#include <cassert>

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
    
};

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

void decompressDfs(unsigned int v, vector <unsigned int> &parentsNewChildren,
                   SuffixTree &compressed, const vector <int> &input) {
    vector <unsigned int> myNewChildren;
    SuffixTree::Node &currentNode = compressed.nodes[v];
    for (auto const &u: currentNode.children) {
        decompressDfs(u, myNewChildren, compressed, input);
    }
    currentNode.children = myNewChildren;
    myNewChildren.clear();
    
    currentNode.lastIndex = min(currentNode.lastIndex * 2, static_cast<unsigned int>(input.size()));
    currentNode.indexOfParentEdge *= 2;
    if (currentNode.leaf != -1) {
        currentNode.leaf *= 2;
    }
    
    if (currentNode.children.size()) {
        vector <unsigned int> similarKids;
        similarKids.push_back(currentNode.children[0]);
        for (unsigned int i = 1; i <= currentNode.children.size(); ++i) {
            if (i != currentNode.children.size()) {
                unsigned int current = currentNode.children[i];
                unsigned int previous = currentNode.children[i - 1];
                while (i < currentNode.children.size() && 
                    input[compressed.nodes[current].indexOfParentEdge] ==
                    input[compressed.nodes[previous].indexOfParentEdge]) {
                    similarKids.push_back(current);
                    ++i;
                    previous = current;
                    if (i < currentNode.children.size()) {
                        current = currentNode.children[i];
                    }
                }
            }
            if (similarKids.size() != 1) {
                unsigned int newNodeIndex = compressed.newNode();
                currentNode = compressed.nodes[v];
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
                    if (currentChild.indexOfParentEdge == currentChild.lastIndex) {
                        assert(currentChild.children.size() == 0);
                        newNode.leaf = currentChild.leaf;
                        compressed.deleteNode(u, newNode.children);
                    }
                }
                myNewChildren.push_back(newNodeIndex);
            } else {
                myNewChildren.push_back(similarKids.back());
            }
            if (i != currentNode.children.size()) {
                similarKids.clear();
                similarKids.push_back(currentNode.children[i]);
            }
        }
    }
    currentNode.children = myNewChildren;
    if (myNewChildren.size() == 1 && v != compressed.root) {
        compressed.deleteNode(v, parentsNewChildren);
    } else {
        parentsNewChildren.push_back(v);
    }
}

void decompress(SuffixTree &compressed, const vector <int> &input) {
    vector <unsigned int> tmp;
    decompressDfs(compressed.root, tmp, compressed, input);
    assert(tmp.size() == 1 && tmp[0] == compressed.root);
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
    consequentPairs = radixSort(radixSort(consequentPairs, [] (const NumberedPair &x) -> int {
        return x.elements.second;
    }), [] (const NumberedPair &x) -> int {
        return x.elements.first;
    });
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
    
    decompress(compressed, input);
    
    return SuffixTree();
}


void dfs_(SuffixTree &tree, int v, int count) {
    for (int i = 0; i < count; ++i) {
        printf("-");
    }
    printf(" %d: %d %d %d %d", v, tree.nodes[v].parent, tree.nodes[v].indexOfParentEdge, tree.nodes[v].lastIndex, tree.nodes[v].leaf);
    for (auto const &u: tree.nodes[v].children) {
        dfs_(tree, u, count + 1);
    }
}

void sample_() {
    SuffixTree sample;
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
    decompress(sample, input);
    dfs_(sample, sample.root, 0);
}


int main() {
    sample_();
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