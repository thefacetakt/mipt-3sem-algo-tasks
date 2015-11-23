#include <cstdio>
#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>
#include <climits>

using std::vector;
using std::make_pair;
using std::max;
using std::min;
using std::pair;

struct SuffixTree {
    struct Node {
        struct Transition {
            unsigned int symbol;
            unsigned int nextIndex;
            Transition(unsigned int symbol, unsigned int nextIndex) : 
                symbol(symbol), 
                nextIndex(nextIndex) {
            }
        };
        int parent;
        unsigned int depth;
        unsigned int indexOfParentEdge;
        int leaf;
        vector <Transition> children;
        
        void addTransition(unsigned int symbol, unsigned int nextIndex) {
            children.push_back(Transition(symbol, newIndex));
        }
    };
    
    vector <Node> nodes;
    vector <unsigned int> pull;
    unsigned int root;
    
    
    SuffixTree() {
        nodes.resize(1);
        root = 0;
        Node &rootNode = nodes[root];
        rootNode.parent = -1;
        rootNode.depth = 0;
        rootNode.indexOfParentEdge = 0;
        rootNode.leaf = 0;
    }
    
    void deleteNode(unsigned int index) {
        Node &node = nodes[index];
        for (auto const &u: node.children) {
            nodes[u].parent = node.parent;
            nodes[u].indexOfParentEdge = node.indexOfParentEdge;
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

void decompressDfs(unsigned int v, vector <SuffixTree::Node::Transition> &parentsNewChildren,
                   SuffixTree &compressed, const vector <int> &input) {
    vector <SuffixTree::Node::Transition> myNewChildren;
    SuffixTree::Node &currentNode = compressed.nodes[v];
    for (auto const &u: currentNode.children) {
        decompressDfs(u.nextIndex, myNewChildren, compressed, length);
    }
    currentNode.children = myNewChildren;
    
    currentNode.depth = min(currentNode.depth * 2, input.size());
    currentNode.indexOfParentEdge *= 2;
    if (currentNode.leaf != -1) {
        currentNode.leaf *= 2;
    }
    if (!currentNode.children.size()) {
        return;
    }
    
    vector <unsigned int> similarKids;
    similarKids.push_back(0);
    for (unsigned int i = 1; i < currentNode.children.size(); ++i) {
        unsigned int current = currentNode.children[i].nextIndex;
        unsigned int previous = currentNode.children[i - 1].nextIndex;
        while (input[compressed.nodes[current].indexOfParentEdge * 2] ==
            input[compressed.nodes[previous].indexOfParentEdge * 2]) {
            ++i;
            unsigned int previous = current;
            unsigned int current = currentNode.children[i].nextIndex;
        }
        unsigned int newNodeIndex = compressed.newNode();
        auto &newNode = compressed.nodes[newNodeIndex];
        similarKids.push_back(i);
    }
    
    
}

void decompress(SuffixTree &compressed, const vector <int> &input) {
    decompressDfs(compressed.root, compressed, input);
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
        newNode.depth = 1;
        newNode.indexOfParentEdge = 0;
        newNode.leaf = 0;
        result.nodes[result.root].addTransition(input[0], newIndex);
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

int main() {
    int n;
    scanf("%d", &n);
    vector <unsigned int> x(n);
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