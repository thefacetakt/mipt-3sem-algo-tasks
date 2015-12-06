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


using std::vector;
using std::make_pair;
using std::max;
using std::min;
using std::pair;
using std::function;
using std::swap;

class SuffixTree {
public:
    class Node {
        vector <unsigned int> children_;
        SuffixTree const *tree;
    public:
        int parent;
        unsigned int indexOfParentEdge;
        unsigned int depth;
        int leaf;

        Node(SuffixTree const *tree=NULL, int parent=-1,
            unsigned int indexOfParentEdge=0, unsigned int depth=0, int leaf=-1
        ) :
            tree(tree),
            parent(parent),
            indexOfParentEdge(indexOfParentEdge),
            depth(depth),
            leaf(leaf) {
        }

        unsigned int lengthOfEdge(unsigned int parentDepth=-1) const {
            if (parentDepth != -1) {
                return depth - parentDepth;
            }
            if (parent != -1) {
                return depth - tree->nodes_[parent].depth;
            }
            return 0;
        }

        unsigned int lastIndex(unsigned int parentDepth=-1) const {
            return indexOfParentEdge + lengthOfEdge(parentDepth);
        }

        unsigned int getFirstHiddenInfo() const {
            return indexOfParentEdge;
        }

        unsigned int getSecondHiddenInfo() const {
            return -(leaf + 2);
        }

        void setHiddenInfo(unsigned int first, unsigned int second) {
            indexOfParentEdge = first;
            leaf = -2 - second;
        }

        unsigned int &operator[](size_t i) {
            return children_[i];
        }

        const unsigned int &operator[](size_t i) const {
            return children_[i];
        }

        void push_back(unsigned int child) {
            children_.push_back(child);
        }

        size_t size() {
            return children_.size();
        }

        void clear() {
            children_.clear();
        }

        vector<unsigned int>::iterator begin() {
            return children_.begin();
        }

        vector<unsigned int>::iterator end() {
            return children_.end();
        }

        vector<unsigned int>::const_iterator begin() const {
            return children_.cbegin();
        }

        vector<unsigned int>::const_iterator end() const {
            return children_.cend();
        }

        void renewChildren(const vector <unsigned int> &newChildren) {
            children_ = newChildren;
        }

        void deleteFirstChild() {
            children_ = vector <unsigned int> (begin() + 1, end());
        }
    };

private:
    vector <Node> nodes_;
    vector <unsigned int> pull_;

public:
    const unsigned int root = 0;

    Node &operator[](size_t i) {
        return nodes_[i];
    }

    const Node &operator[](size_t i) const {
        return nodes_[i];
    }

    SuffixTree() {
        nodes_.resize(1);
    }

    inline void checkNode(const Node &node) const {
        #ifdef _DEBUG
            assert(node.parent == -1 || nodes_[node.parent].depth < node.depth);
        #endif
    }

    void deleteUselessNode(unsigned int v, unsigned int inParentIndex,
        int leaf
     ) {
        unsigned int parent = nodes_[v].parent;
        nodes_[parent][inParentIndex] = -1;

        #ifdef _DEBUG
            assert(nodes_[v].size() <= 1);
        #endif

        for (auto const &u: nodes_[v]) {
            nodes_[parent][inParentIndex] = u;
            nodes_[u].parent = parent;
            nodes_[u].indexOfParentEdge = leaf + nodes_[parent].depth;
            checkNode(nodes_[u]);
        }
        nodes_[v].clear();
        pull_.push_back(v);
    }

    unsigned int newNode(
        int parent=-1,
        unsigned int indexOfParentEdge=0,
        unsigned int depth=0,
        int leaf=-1
    ) {
        Node resultNode = Node(this, parent, indexOfParentEdge, depth, leaf);
        if (pull_.size()) {
            unsigned int returnValue = pull_.back();
            pull_.pop_back();
            nodes_[returnValue] = resultNode;
            return returnValue;
        }
        nodes_.push_back(resultNode);
        return nodes_.size() - 1;
    }

    unsigned int splitEdge(unsigned int parent, unsigned int childIndex,
        unsigned int length
    ) {
        unsigned int child = nodes_[parent][childIndex];
        unsigned int newNodeIndex = newNode(parent,
            nodes_[child].indexOfParentEdge, nodes_[parent].depth + length);
        Node &newNode = nodes_[newNodeIndex];
        nodes_[child].indexOfParentEdge = newNode.lastIndex();
        nodes_[child].parent = newNodeIndex;
        newNode.push_back(child);
        nodes_[parent][childIndex] = newNodeIndex;

        checkNode(nodes_[newNodeIndex]);

        return newNodeIndex;
    }
};


void buildEulerDfs(unsigned int v, unsigned int &count, unsigned int depth,
    const SuffixTree &tree, vector <EulerPair> euler,
    vector <unsigned int> realNode
) {
    auto const &currentNode = tree[v];
    realNode.push_back(v);
    euler.push_back(EulerPair(depth, count));
    int currentCount = count;
    for (auto const &u: currentNode) {
        buildEulerDfs(u, ++count, depth + 1, tree, euler, realNode);
        euler.push_back(EulerPair(depth, count));
    }
}

void buildEuler(const SuffixTree &tree, vector <EulerPair> euler,
    vector <unsigned int> realNode
) {
    unsigned int count = 0;
    buildEulerDfs(tree.root, count, 0, tree, euler, realNode);
}

class RandomLCPGetter {
private:
    const SuffixTree &tree;
    vector <EulerPair> euler;
    vector <unsigned int> realNode;
    vector <unsigned int> irrealNode;
    LCA lcaGetter;
public:
    RandomLCPGetter(const SuffixTree &tree, unsigned int inputLength) :
        tree(tree)
    {
        buildEuler(tree, euler, realNode);
        irrealNode.resize(inputLength + 1);
        for (unsigned int i = 0; i < realNode.size(); ++i) {
            int leaf = tree[realNode[i]].leaf;
            if (leaf != -1) {
                irrealNode[leaf] = i;
            }
        }
        lcaGetter = LCA(euler);
    }

    unsigned int lca(unsigned int i, unsigned int j) const {
        return realNode[lcaGetter.lca(irrealNode[i], irrealNode[j])];
    }
    unsigned int lcp(unsigned int i, unsigned int j) const {
        return tree[lca(i, j)].depth;
    }
};




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

vector <int> compressInput(const vector <int> &input) {
    vector <NumberedPair> consequentPairs;
    for (unsigned int i = 0; i < input.size(); i += 2) {
        consequentPairs.push_back(NumberedPair(
                input[i],
                (i + 1 == input.size() ? -1: input[i + 1]),
                i / 2
            )
        );
    }
    consequentPairs = radixSort(
        radixSort(consequentPairs, secondElementOfNumberedPair),
        firstElementOfNumberedPair
    );

    vector <int> pairedString((input.size() + 1) / 2);
    unsigned int currentNumber = 0;

    for (unsigned int i = 0; i < consequentPairs.size(); ++i) {
        if (i && consequentPairs[i].elements
                != consequentPairs[i - 1].elements) {
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

    return pairedString;
}

int decompressDfs(unsigned int v, unsigned int inParentIndex, SuffixTree &tree,
    const vector <int> &input, unsigned int depth
) {
    int leaf = -1;
    unsigned int oldDepth = tree[v].depth;
    tree[v].depth = depth;

    for (unsigned int i = 0; i < tree[v].size(); ++i) {
        unsigned int u = tree[v][i];
        int newLeaf = decompressDfs(u, i, tree, input,
            depth + tree[u].lengthOfEdge(oldDepth) * 2
            - (tree[v].lastIndex(oldDepth) == input.size() / 2 + 1)
        );
        if (leaf == -1) {
            leaf = newLeaf;
        }
    }
    tree[v].indexOfParentEdge *= 2;
    if (tree[v].leaf != -1) {
        tree[v].leaf = min(tree[v].leaf * 2, static_cast<int>(input.size()));
    }

    vector <unsigned int> myNewChildren;
    if (tree[v].size()) {
        vector <unsigned int> similarKids;
        unsigned int firstKidIndex = 0;
        similarKids.push_back(tree[v][0]);
        for (unsigned int i = 1; i <= tree[v].size(); ++i) {
            if (i != tree[v].size()) {
                unsigned int current = tree[v][i];
                unsigned int previous = tree[v][i - 1];
                while (i < tree[v].size()
                    && input[tree[current].indexOfParentEdge]
                        == input[tree[previous].indexOfParentEdge]
                ) {
                    similarKids.push_back(current);
                    ++i;
                    previous = current;
                    if (i != tree[v].size()) {
                        current = tree[v][i];
                    }
                }
            }
            if (similarKids.size() != 1) {
                unsigned int newNodeIndex = tree.splitEdge(v, firstKidIndex, 1);
                auto &newNode = tree[newNodeIndex];
                for (unsigned int j = 1; j < similarKids.size(); ++j) {
                    auto &currentChild = tree[similarKids[j]];
                    currentChild.parent = newNodeIndex;
                    ++currentChild.indexOfParentEdge;
                    newNode.push_back(similarKids[j]);
                }
                if (newNode.size() && tree[newNode[0]].lengthOfEdge() == 0) {
                    newNode.leaf = tree[newNode[0]].leaf;
                    tree.deleteUselessNode(newNode[0], 0, newNode.leaf);
                    newNode.deleteFirstChild();
                }
                myNewChildren.push_back(newNodeIndex);

                tree.checkNode(tree[newNodeIndex]);
            } else {
                myNewChildren.push_back(similarKids.back());
            }
            if (i != tree[v].size()) {
                similarKids.clear();
                similarKids.push_back(tree[v][i]);
                firstKidIndex = i;
            }
        }
    }
    tree[v].renewChildren(myNewChildren);
    if (tree[v].size() == 1 && v != tree.root && tree[v].leaf == -1) {
        tree.deleteUselessNode(v, inParentIndex, leaf);
    }
    if (tree[v].leaf != -1) {
        leaf = tree[v].leaf;
    }
    return leaf;
}

void decompress(SuffixTree &tree, const vector <int> &input) {
    decompressDfs(tree.root, -1, tree, input, 0);
}

void checkDfs(const SuffixTree &tree, unsigned int v, const vector<int> &input
) {
    #ifdef _PRINT_DBG
        printf("%d -> %d [label=\"", tree[v].parent, v);
        for (int i = tree[v].indexOfParentEdge;
            i != tree[v].lastIndex(); ++i) {
            printf("%d", input[i]);
        }
        printf(", %d, %d\"]\n", tree[v].leaf, tree[v].depth);
    #endif

    tree.checkNode(tree[v]);
    for (auto const &u: tree[v]) {
        checkDfs(tree, u, input);
    }
}

void printHiddenDfs(const SuffixTree &tree, unsigned int v) {
    printf("%d -> %d [label=\"%d:%d, %d %d\"]\n", tree[v].parent, v,
        tree[v].indexOfParentEdge, tree[v].leaf, tree[v].depth
    );
    for (auto const &u: tree[v]) {
        printHiddenDfs(tree, u);
    }
}

void suffixArrayDfs(unsigned int v, const SuffixTree &tree,
    vector <unsigned int> &suffix
) {
    auto &currentNode = tree[v];
    if (currentNode.leaf != -1) {
        suffix.push_back(currentNode.leaf);
    }

    for (auto const &u: currentNode) {
        suffixArrayDfs(u, tree, suffix);
    }
}

void buildSuffixArray(const SuffixTree &tree, vector <unsigned int> &suffix) {
    suffixArrayDfs(tree.root, tree, suffix);
}

vector <unsigned int> buildOddSuffixArray(
    const vector <unsigned int> &evenSuffix,
    const vector <int> &input
) {
    vector <unsigned int> antiSuffixArray(input.size());
    for (unsigned int i = 0; i < evenSuffix.size(); ++i) {
        antiSuffixArray[evenSuffix[i]] = i;
    }

    vector <unsigned int> oddSuffix;
    for (unsigned int i = 1; i < input.size(); i += 2) {
        oddSuffix.push_back(i);
    }

    oddSuffix = radixSort(radixSort(oddSuffix, [&input, &antiSuffixArray]
        (unsigned int i) -> int {
            return (i + 1 == input.size() ? -1 : antiSuffixArray[i + 1]);
        }),
    [&input] (unsigned int i) -> unsigned int {
        return input[i];
    });

    return oddSuffix;
}

vector <unsigned int> buildOddLcp(const RandomLCPGetter &getter,
    const vector <unsigned int> &oddSuffix,
    const vector <int> &input
) {
    vector <unsigned int> oddLcp(oddSuffix.size() - 1);
    for (unsigned int i = 0; i + 1 < oddSuffix.size(); ++i) {
        if (input[oddSuffix[i]] == input[oddSuffix[i + 1]]) {
            if (oddSuffix[i + 1] + 1 == input.size()
                || oddSuffix[i] + 1 == input.size()) {
                oddLcp[i] = 1;
            } else {
                oddLcp[i]
                    = getter.lcp(oddSuffix[i] + 1, oddSuffix[i + 1] + 1) + 1;
            }
        } else {
            oddLcp[i] = 0;
        }
    }

    #ifdef _DEBUG
        for (unsigned int i = 0; i < oddSuffix.size(); ++i) {
            if (i + 1 != oddSuffix.size()) {
                assert(
                    vector<int>(input.begin() + oddSuffix[i],
                        input.begin() + oddSuffix[i] + oddLcp[i])
                    ==
                    vector<int>(input.begin() + oddSuffix[i + 1],
                        input.begin() + oddSuffix[i + 1] + oddLcp[i])
                    &&
                        (oddSuffix[i] + oddLcp[i] == input.size()
                        ||
                        oddSuffix[i + 1] + oddLcp[i] == input.size()
                        ||
                    vector<int>(input.begin() + oddSuffix[i],
                        input.begin() + oddSuffix[i] + oddLcp[i] + 1)
                    !=
                    vector<int>(input.begin() + oddSuffix[i + 1],
                        input.begin() + oddSuffix[i + 1] + oddLcp[i] + 1))
                );
            }
            #ifdef _PRINT_DBG
                for (unsigned int j = oddSuffix[i]; j < input.size(); ++j) {
                    printf("%d ", input[j]);
                }
                printf("| %d\n", (i == oddLcp.size() ? -1 : oddLcp[i]));
            #endif
        }
    #endif

    return oddLcp;
}

SuffixTree buildSuffixTreeFromSA(vector <unsigned int> &sa,
    vector <unsigned int> &lcp,
    unsigned int length
) {
    vector <int> tmp;
    SuffixTree result = buildTempSuffixTree(tmp);
    int newNodeIndex
        = result.newNode(result.root, sa[0], length - sa[0], sa[0]);
    result[result.root].push_back(newNodeIndex);
    unsigned int current = newNodeIndex;

    for (unsigned int i = 1; i < sa.size(); ++i) {
        while (result[current].parent != -1
            && result[result[current].parent].depth >= lcp[i - 1]
        ) {
            current = result[current].parent;
        }
        unsigned int parent;
        if (result[current].parent != -1 &&
            result[result[current].parent].depth == lcp[i - 1]
        ) {
            parent = result[current].parent;
        } else if (result[current].depth == lcp[i - 1]) {
            parent = current;
        } else {
            unsigned int currentParent = result[current].parent;
            parent = result.splitEdge(currentParent,
                result[currentParent].size() - 1,
                lcp[i - 1] - result[currentParent].depth
            );
            result[parent].leaf = (length - sa[i] == lcp[i - 1] ? sa[i] : -1);
        }
        if (lcp[i - 1] != length - sa[i]) {
            newNodeIndex = result.newNode(parent, sa[i] + lcp[i - 1],
                length - sa[i], sa[i]
            );
            result[parent].push_back(newNodeIndex);
            parent = newNodeIndex;
        }
        current = parent;
    }
    return result;
}

SuffixTree buildOddSuffixTree(const SuffixTree &even,
    const vector<int> &input
) {
    vector <unsigned int> evenSuffix;
    buildSuffixArray(even, evenSuffix);
    RandomLCPGetter evenGetter(even, input.size());

    vector <unsigned int> oddSuffix = buildOddSuffixArray(evenSuffix, input);
    vector <unsigned int> oddLcp = buildOddLcp(evenGetter, oddSuffix, input);
    return buildSuffixTreeFromSA(oddSuffix, oddLcp, input.size());
}

inline void copyNodeExceptParentAndChildren(const SuffixTree::Node &from,
    SuffixTree::Node &to
) {
    to.depth = from.depth;
    to.indexOfParentEdge = from.indexOfParentEdge;
    to.leaf = from.leaf;
}

void copySubTree(const SuffixTree &from, SuffixTree &to,
    unsigned int fromStart, unsigned int toStart
);

unsigned int appendCopyNode(const SuffixTree &from, SuffixTree &to,
    unsigned int toStart, unsigned int u
) {
    unsigned int newStart = to.newNode(toStart);
    to[newStart].parent = toStart;
    to[toStart].push_back(newStart);
    copyNodeExceptParentAndChildren(from[u], to[newStart]);
    copySubTree(from, to, u, newStart);
    return newStart;
}

void copySubTree(const SuffixTree &from, SuffixTree &to,
    unsigned int fromStart, unsigned int toStart
) {
    for (auto const &u: from[fromStart]) {
        appendCopyNode(from, to, toStart, u);
    }
}


struct MergeNodesStruct {
    SuffixTree &tree;
    const vector <int> &input;
    unsigned int v;
    unsigned int childIndex;
    unsigned int child;
    unsigned int length;
    int symbol;

    MergeNodesStruct(SuffixTree &tree,
        unsigned int v, const vector <int> &input
    ) : tree(tree), v(v), input(input) {
        childIndex = 0;
    }

    bool end() const {
        return child == tree[v].size();
    }

    void evaluate() {
        child = tree[v][childIndex];
        symbol = input[tree[child].indexOfParentEdge];
    }

    void split(unsigned int splitLength) {
        child = tree.splitEdge(v, childIndex, splitLength);
    }
};

template<class ActionType>
void doSomething(MergeNodesStruct merging[2], ActionType action) {
    action(merging[0]);
    action(merging[1]);
}

template<class Compare>
MergeNodesStruct &minimal(MergeNodesStruct merging[2], Compare comp) {
    if (comp(merging[0]) < comp(merging[1])) {
        return merging[0];
    }
    return merging[1];
}

void mergeNodes(unsigned int first, unsigned int second, unsigned int to,
    SuffixTree &result, SuffixTree &tree1, SuffixTree &tree2,
    const vector <int> &input
) {
    MergeNodesStruct merging[] = {MergeNodesStruct(tree1, first, input),
        MergeNodesStruct(tree2, second, input)
    };
    while (!merging[0].end() && !merging[1].end()) {
        doSomething(merging, [] (MergeNodesStruct &instance) {
            instance.evaluate();
        });
        if (merging[0].symbol == merging[1].symbol) {
            doSomething(merging, [] (MergeNodesStruct &instance) {
                instance.length = instance.tree[instance.child].lengthOfEdge();
            });
            MergeNodesStruct &splitVictim = minimal(merging,
                [] (const MergeNodesStruct &instance) -> int {
                    return -static_cast<int>(instance.length);
                }
            );
            unsigned int minimalLength
                = min(merging[0].length, merging[1].length);
            splitVictim.split(minimalLength);

            unsigned int newNodeIndex = result.newNode(to);
            auto &newNode = result[newNodeIndex];
            result[to].push_back(newNodeIndex);

            if (minimalLength == 1) {
                MergeNodesStruct &recipient = minimal(merging,
                    [] (const MergeNodesStruct &instance)
                        -> pair<unsigned int, int> {
                        return make_pair(instance.length,
                            -instance.tree[instance.child].leaf);
                    }
                );
                copyNodeExceptParentAndChildren(recipient.tree[recipient.child],
                    newNode
                );
            } else {
                newNode.setHiddenInfo(merging[0].child, merging[1].child);
                newNode.depth = minimalLength + result[to].depth;
            }
            mergeNodes(merging[0].child, merging[1].child, newNodeIndex,
                result, tree1, tree2, input);
            doSomething(merging, [] (MergeNodesStruct &instance) {
                ++instance.childIndex;
            });
        } else {
            MergeNodesStruct &recipient = minimal(merging,
                [] (const MergeNodesStruct &instance) -> int {
                    return instance.symbol;
                }
            );
            ++recipient.childIndex;
            appendCopyNode(recipient.tree, result, to, recipient.child);
        }
    }
    MergeNodesStruct &recipient = minimal(merging,
        [] (const MergeNodesStruct & instance) -> bool {
            return instance.end();
        }
    );
    for (; !recipient.end(); ++recipient.childIndex) {
        appendCopyNode(recipient.tree, result, to, recipient.child);
    }
}

bool findSuffixesDfs(unsigned int v, const SuffixTree &tree,
    vector <unsigned int> &output
) {
    while (output.size() <= v) {
        output.push_back(-1);
    }
    output[v] = tree[v].leaf;
    for (auto const &u: tree[v]) {
        if (findSuffixesDfs(u, tree, output)) {
            if (output[v] == -1) {
                output[v] = output[u];
            }
        }
    }
    return (output[v] != -1);
}

vector <unsigned int> findSuffixes(const SuffixTree &tree) {
    vector <unsigned int> result;
    findSuffixesDfs(tree.root, tree, result);
    return result;
}

SuffixTree mergeTrees(SuffixTree &tree1, SuffixTree &tree2,
    const vector <int> &input
) {
    vector <int> tmp;
    SuffixTree merged = buildTempSuffixTree(tmp);
    merged[merged.root].leaf = -1;
    mergeNodes(tree1.root, tree2.root, merged.root, merged,
        tree1, tree2, input
    );

    #ifdef _PRINT_DBG
        printf("MERGED\n");
        printHiddenDfs(merged, merged.root);
    #endif

    #ifdef _DEBUG
        #ifdef _PRINT_DBG
            printf("NEW TREE1\n");
        #endif
        checkDfs(tree1, tree1.root, input);
        #ifdef _PRINT_DBG
            printf("NEW TREE2\n");
        #endif
        checkDfs(tree2, tree2.root, input);
    #endif


    vector <unsigned int> firstSuffixes = findSuffixes(tree1);
    vector <unsigned int> secondSuffixes = findSuffixes(tree2);

    vector <unsigned int> trueLength = computeTrueLength(merged, tree1, tree2, input.size());
    correctMerge(merged.root, -1, merged, tree1, tree2, trueLength, input, firstSuffixes, secondSuffixes);

    return merged;
}

SuffixTree buildTempSuffixTree(const vector <int> &input) {
    if (input.size() == 0) {
        return SuffixTree();
    }
    if (input.size() == 1) {
        SuffixTree result;
        unsigned int newNodeIndex = result.newNode(result.root, 0, 1, 0);
        result[result.root].push_back(newNodeIndex);
        return result;
    }

    SuffixTree compressed = buildTempSuffixTree(compressInput(input));
    decompress(compressed, input);
    SuffixTree &even = compressed;

    #ifdef _DEBUG
        #ifdef _PRINT_DBG
            printf("EVEN:\n");
        #endif
        checkDfs(even, even.root, input);
    #endif

    SuffixTree odd = buildOddSuffixTree(even, input);

    #ifdef _DEBUG
        #ifdef _PRINT_DBG
            printf("ODD:\n");
        #endif
        checkDfs(odd, odd.root, input);
    #endif

    SuffixTree almostResult = mergeTrees(even, odd, input);
}
