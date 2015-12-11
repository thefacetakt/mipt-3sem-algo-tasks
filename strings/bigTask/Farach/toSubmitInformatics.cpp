#ifndef _MINIMAL_PAIR
#define _MINIMAL_PAIR

template<typename T>
struct MinimalPair {
    T element;
    unsigned int index;

    MinimalPair() {
    }

    MinimalPair(T element, unsigned int index) :
        element(element), index(index) {
    }

    bool operator<(const MinimalPair &other) const {
        return (element == other.element ?
            index < other.index : element < other.element
        );
    }
};

#endif

#ifndef _SPARSE_TABLE
#define _SPARSE_TABLE


#include <vector>
#include <algorithm>
#include <climits>

using std::vector;
using std::min;

struct SparseTable {
private:
    vector<vector<MinimalPair<unsigned int> > > st;
    vector<unsigned int> fastLog;
    unsigned int n;

    void init();

public:
    SparseTable(const vector<unsigned int> &elements);
    SparseTable();

    unsigned int minimum(unsigned int i, unsigned int j) const;

    MinimalPair<unsigned int> operator[](unsigned int i) const;
};

#endif

#include <vector>
#include <algorithm>
#include <climits>

using std::vector;
using std::min;


SparseTable::SparseTable() {
}

SparseTable::SparseTable(const vector<unsigned int> &elements) {
    n = elements.size();
    fastLog.resize(n + 1);
    fastLog[1] = 0;
    for (unsigned int i = 2; i <= n; ++i) {
        fastLog[i] = fastLog[i / 2] + 1;
    }
    st.resize(fastLog[n] + 1);
    st[0].resize(n);
    for (unsigned int i = 0; i < n; ++i) {
        st[0][i] = MinimalPair<unsigned int>(elements[i], i);
    }
    init();
}

void SparseTable::init() {
    for (unsigned int k = 1; k <= fastLog[n]; ++k) {
        st[k].resize(n - (1 << k) + 1);
        for (unsigned int i = 0; i < st[k].size(); ++i) {
            st[k][i] = min(st[k - 1][i], st[k - 1][i + (1 << (k - 1))]);
        }
    }
}

unsigned int SparseTable::minimum(unsigned int i, unsigned int j) const {
    if (j < i) {
        return UINT_MAX;
    }
    unsigned int length = (j - i + 1);
    return min(st[fastLog[length]][i],
        st[fastLog[length]][j - (1 << fastLog[length]) + 1]).index;
}

MinimalPair<unsigned int> SparseTable::operator[](unsigned int i) const {
   return st[0][i];
}

#ifndef _RMQpm1
#define _RMQpm1

#include <vector>
#include <algorithm>
#include <climits>
#include <cstdio>

using std::vector;
using std::max;
using std::min;

struct RMQpm1 {
private:
    SparseTable st;
    vector <unsigned int> elements;
    vector <vector <MinimalPair<unsigned int> > > prefixMins;
    vector <vector <MinimalPair<unsigned int> > > suffixMins;
    unsigned int n;
    unsigned int block;
    vector <vector<vector<MinimalPair<int> > > > dp;
    vector <unsigned int> type;

public:
    RMQpm1();
    RMQpm1(const vector<unsigned int> &elements);

    unsigned int minimum(unsigned int i, unsigned int j) const;
};

#endif

#include <vector>
#include <algorithm>
#include <climits>
#include <cstdio>

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

#ifndef _EULER_PAIR
#define _EULER_PAIR

struct EulerPair {
    unsigned int depth;
    unsigned int node;

    EulerPair();

    EulerPair(unsigned int depth, unsigned int node);
};

#endif


EulerPair::EulerPair() {
}

EulerPair::EulerPair(unsigned int depth, unsigned int node) :
    depth(depth),
    node(node) {
}

#ifndef _LCA
#define _LCA

#include <vector>
#include <algorithm>
#include <utility>
#include <climits>
#include <cstdio>


using std::vector;
using std::max;
using std::min;

class LCA {
private:
    vector <unsigned int> myEuler;
    vector <unsigned int> first;
    vector <unsigned int> last;

    RMQpm1 rmq;

public:

    LCA();

    LCA(const vector<EulerPair> &euler);

    unsigned int lca(unsigned int u, unsigned int v) const;
};

#endif


#ifndef _LCA_CPP
#define _LCA_CPP

#include <vector>
#include <algorithm>
#include <utility>
#include <climits>
#include <cstdio>



using std::vector;
using std::max;
using std::min;

LCA::LCA() {
}

LCA::LCA(const vector<EulerPair> &euler) {
    unsigned int n = (euler.size() + 1) / 2;
    vector <unsigned int> toRMQ(euler.size());
    myEuler.resize(euler.size());
    first.assign(n, -1);
    last.assign(n, -1);
    for (unsigned int i = 0; i < euler.size(); ++i) {
        if (first[euler[i].node] == -1) {
            first[euler[i].node] = i;
        }
        last[euler[i].node] = i;
        toRMQ[i] = euler[i].depth;
        myEuler[i] = euler[i].node;
    }

    rmq = RMQpm1(toRMQ);
}


unsigned int LCA::lca(unsigned int u, unsigned int v) const {
    return myEuler[rmq.minimum(min(first[u], first[v]), max(last[u], last[v]))];
}

#endif

#ifndef _SUFFIX_TREE
#define _SUFFIX_TREE

#include <vector>
#include <cassert>

using std::vector;

class SuffixTree {
public:
    class Node {
        vector <unsigned int> children_;
    public:
        int parent;
        unsigned int indexOfParentEdge;
        unsigned int depth;
        int leaf;

        Node(int parent=-1,
            unsigned int indexOfParentEdge=0, unsigned int depth=0, int leaf=-1
        );
        unsigned int lengthOfEdge(const SuffixTree &tree,
            unsigned int parentDepth=-1
        ) const;

        unsigned int lastIndex(const SuffixTree &tree,
            unsigned int parentDepth=-1
        ) const;

        unsigned int getFirstHiddenInfo() const;

        unsigned int getSecondHiddenInfo() const;

        unsigned int getHiddenInfo(unsigned int i) const;

        void setHiddenInfo(unsigned int first, unsigned int second);

        unsigned int &operator[](size_t i);

        const unsigned int &operator[](size_t i) const;

        void push_back(unsigned int child);

        size_t size() const;

        void clear();

        vector<unsigned int>::iterator begin();

        vector<unsigned int>::iterator end();

        vector<unsigned int>::const_iterator begin() const;

        vector<unsigned int>::const_iterator end() const;

        void renewChildren(const vector <unsigned int> &newChildren);

        void deleteFirstChild();

        void resize(size_t size);

        bool isHiddenInfo() const;
    };

private:
    vector <Node> nodes_;
    vector <unsigned int> pull_;

public:
    const unsigned int root = 0;

    Node &operator[](size_t i);

    const Node &operator[](size_t i) const;

    SuffixTree();

    void checkNode(const Node &node) const;

    void deleteUselessNode(unsigned int v, unsigned int inParentIndex,
        int leaf
    );

    unsigned int newNode(
        int parent=-1,
        unsigned int indexOfParentEdge=0,
        unsigned int depth=0,
        int leaf=-1
    );

    unsigned int splitEdge(unsigned int parent, unsigned int childIndex,
        unsigned int length
    );
};

#endif

#include <vector>
#include <cassert>

using std::vector;

SuffixTree::Node::Node(int parent,
    unsigned int indexOfParentEdge, unsigned int depth, int leaf
) :
    parent(parent),
    indexOfParentEdge(indexOfParentEdge),
    depth(depth),
    leaf(leaf) {
}

unsigned int SuffixTree::Node::lengthOfEdge(const SuffixTree &tree,
    unsigned int parentDepth
) const {
    if (parentDepth != -1) {
        return depth - parentDepth;
    }
    if (parent != -1) {
        return depth - tree[parent].depth;
    }
    return depth;
}

unsigned int SuffixTree::Node::lastIndex(const SuffixTree &tree,
    unsigned int parentDepth
) const {
    return indexOfParentEdge + lengthOfEdge(tree, parentDepth);
}

unsigned int SuffixTree::Node::getFirstHiddenInfo() const {
    return indexOfParentEdge;
}

unsigned int SuffixTree::Node::getSecondHiddenInfo() const {
    return -(leaf + 2);
}

unsigned int SuffixTree::Node::getHiddenInfo(unsigned int i) const {
    #ifdef _GLIBCXX_DEBUG
        assert(0 <= i && i <= 1);
    #endif
    return (i == 0 ? getFirstHiddenInfo() : getSecondHiddenInfo());
}

void SuffixTree::Node::setHiddenInfo(unsigned int first, unsigned int second) {
    indexOfParentEdge = first;
    leaf = -2 - second;
}

unsigned int &SuffixTree::Node::operator[](size_t i) {
    return children_[i];
}

const unsigned int &SuffixTree::Node::operator[](size_t i) const {
    return children_[i];
}

void SuffixTree::Node::push_back(unsigned int child) {
    children_.push_back(child);
}

size_t SuffixTree::Node::size() const {
    return children_.size();
}

void SuffixTree::Node::clear() {
    children_.clear();
}

vector<unsigned int>::iterator SuffixTree::Node::begin() {
    return children_.begin();
}

vector<unsigned int>::iterator SuffixTree::Node::end() {
    return children_.end();
}

vector<unsigned int>::const_iterator SuffixTree::Node::begin() const {
    return children_.cbegin();
}

vector<unsigned int>::const_iterator SuffixTree::Node::end() const {
    return children_.cend();
}

void SuffixTree::Node::renewChildren(const vector <unsigned int> &newChildren) {
    children_ = newChildren;
}

void SuffixTree::Node::deleteFirstChild() {
    children_ = vector <unsigned int> (begin() + 1, end());
}

void SuffixTree::Node::resize(size_t size) {
    children_.resize(size);
}

bool SuffixTree::Node::isHiddenInfo() const {
    return leaf <= -2;
}

SuffixTree::Node &SuffixTree::operator[](size_t i) {
    return nodes_[i];
}

const SuffixTree::Node &SuffixTree::operator[](size_t i) const {
    return nodes_[i];
}

SuffixTree::SuffixTree() {
    nodes_.resize(1);
}

void SuffixTree::checkNode(const Node &node) const {
    #ifdef _DEBUG
        assert(node.parent == -1 || nodes_[node.parent].depth < node.depth);
    #endif
}

void SuffixTree::deleteUselessNode(unsigned int v, unsigned int inParentIndex,
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

unsigned int SuffixTree::newNode(
    int parent,
    unsigned int indexOfParentEdge,
    unsigned int depth,
    int leaf
) {
    Node resultNode = Node(parent, indexOfParentEdge, depth, leaf);
    if (pull_.size()) {
        unsigned int returnValue = pull_.back();
        pull_.pop_back();
        nodes_[returnValue] = resultNode;
        return returnValue;
    }
    nodes_.push_back(resultNode);
    return nodes_.size() - 1;
}

unsigned int SuffixTree::splitEdge(unsigned int parent, unsigned int childIndex,
    unsigned int length
) {
    unsigned int child = nodes_[parent][childIndex];
    unsigned int newNodeIndex = newNode(parent,
        nodes_[child].indexOfParentEdge, nodes_[parent].depth + length);
    Node &newNode = nodes_[newNodeIndex];
    nodes_[child].indexOfParentEdge = newNode.lastIndex(*this);
    nodes_[child].parent = newNodeIndex;
    newNode.push_back(child);
    nodes_[parent][childIndex] = newNodeIndex;

    checkNode(nodes_[newNodeIndex]);

    return newNodeIndex;
}

class RandomLCPGetter;
struct NumberedPair;

template<class T>
class IndexedPair;

struct MergeNodesStruct;
struct MergeTreesStruct;

#ifndef _FARACH
#define _FARACH

#include <cstdio>
#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>
#include <climits>
#include <cassert>
#include <functional>
#include <string>

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


#include <cstdio>
#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>
#include <climits>
#include <cassert>
#include <functional>
#include <string>


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
) {
    auto const &currentNode = tree[v];
    realNode.push_back(v);
    euler.push_back(EulerPair(depth, count));
    int currentCount = count;
    for (auto const &u: currentNode) {
        buildEulerDfs(u, ++count, depth + 1, tree, euler, realNode);
        euler.push_back(EulerPair(depth, currentCount));
    }
}

void buildEuler(const SuffixTree &tree, vector <EulerPair> &euler,
    vector <unsigned int> &realNode
) {
    unsigned int count = 0;
    buildEulerDfs(tree.root, count, 0, tree, euler, realNode);
}

RandomLCPGetter usualGetter(const SuffixTree &tree, unsigned int inputLength) {
    return RandomLCPGetter(tree, inputLength,
        [&tree](ELeafGetterMode mode,
            unsigned int vertice=-1, unsigned int number=-1
        ) -> int {
            if (mode == NUMBER) {
                return 1;
            }
            return tree[vertice].leaf;
        }
    );
}

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
            depth + tree[u].lengthOfEdge(tree, oldDepth) * 2
            - (tree[u].lastIndex(tree, oldDepth) == input.size() / 2 + 1)
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
                if (newNode.size()
                    && tree[newNode[0]].lengthOfEdge(tree) == 0
                ) {
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
            i != tree[v].lastIndex(tree); ++i) {
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
        tree[v].indexOfParentEdge, tree[v].lastIndex(tree),
        tree[v].leaf, tree[v].depth
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
    RandomLCPGetter evenGetter = usualGetter(even, input.size());

    vector <unsigned int> oddSuffix = buildOddSuffixArray(evenSuffix, input);
    vector <unsigned int> oddLcp = buildOddLcp(evenGetter, oddSuffix, input);
    return buildSuffixTreeFromSA(oddSuffix, oddLcp, input.size());
}

void copyNodeExceptParentAndChildren(const SuffixTree::Node &from,
    SuffixTree::Node &to
) {
    to.depth = from.depth;
    to.indexOfParentEdge = from.indexOfParentEdge;
    to.leaf = from.leaf;
}

unsigned int appendCopyNode(const SuffixTree &from, SuffixTree &to,
    unsigned int toStart, unsigned int u
) {
    unsigned int newStart = to.newNode(toStart);
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

void mergeNodes(unsigned int first, unsigned int second, unsigned int to,
    SuffixTree &result, SuffixTree &tree1, SuffixTree &tree2,
    const vector <int> &input
) {
    IndexedPair<MergeNodesStruct> merging(MergeNodesStruct(tree1, first, input),
        MergeNodesStruct(tree2, second, input)
    );
    while (!merging[0].end() && !merging[1].end()) {
        doSomething(merging, [] (MergeNodesStruct &instance) {
            instance.evaluate();
        });
        if (merging[0].symbol == merging[1].symbol) {
            doSomething(merging, [] (MergeNodesStruct &instance) {
                instance.length
                    = instance.tree[instance.child].lengthOfEdge(instance.tree);
            });
            unsigned int minimalLength
                = min(merging[0].length, merging[1].length);
            if (merging[0].length != merging[1].length) {
                MergeNodesStruct &splitVictim = minimal(merging,
                    [] (const MergeNodesStruct &instance) -> int {
                        return -static_cast<int>(instance.length);
                    }
                );

                splitVictim.split(minimalLength);
            }
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
        recipient.evaluate();
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

void findLcaSuffix(unsigned int v, const SuffixTree &merged,
    IndexedPair<const SuffixTree &> trees,
    vector <IndexedPair<unsigned int> > &output
) {
    while (output.size() <= v) {
        output.push_back(IndexedPair<unsigned int>(-1, -1));
    }
    if (merged[v].leaf >= 0) {
        output[v][merged[v].leaf % 2] = merged[v].leaf;
    } else if (merged[v].isHiddenInfo()) {
        for (unsigned int i = 0; i < 2; ++i) {
            output[v][i] = trees[i][merged[v].getHiddenInfo(i)].leaf;
        }
    }
    int reserve = -1;
    for (auto const &u: merged[v]) {
        findLcaSuffix(u, merged, trees, output);
        if (output[u][0] != -1) {
            if (output[v][0] == -1) {
                output[v][0] = output[u][0];
                if (output[u][1] != -1) {
                    reserve = output[u][1];
                }
            } else if (reserve != -1) {
                output[v][0] = output[u][0];
                output[v][1] = reserve;
            }
        }
        if (output[u][1] != -1 && output[u][1] != reserve) {
            if (output[v][1] == -1) {
                output[v][1] = output[u][1];
            }
        }
    }
    #ifdef _DEBUG
        assert(!merged[v].isHiddenInfo()
            || (output[v][0] != -1 && output[v][1] != -1));
    #endif
}

void trueLengthFillDfs(unsigned int v, const SuffixTree &merged,
    const RandomLCPGetter &getter, IndexedPair<unsigned int> suffix,
    unsigned int length,
    vector <unsigned int> &output
) {
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
    ++suffix[0], ++suffix[1];
    if (suffix[0] == length || suffix[1] == length) {
        output[v] = 1;
        return;
    }
    unsigned int u = getter.lca(suffix[0], suffix[1]);
    trueLengthFillDfs(u, merged, getter, suffix, length, output);
    output[v] = output[u] + 1;
    #ifdef _DEBUG
        assert(output[v] <= merged[v].depth);
    #endif

}

vector <unsigned int> computeTrueLength(const SuffixTree &merged,
    IndexedPair<const SuffixTree &> trees, unsigned int inputLength
) {
    RandomLCPGetter getter(merged, inputLength,
        [&merged, &trees] (ELeafGetterMode mode,
            unsigned int vertice=-1, unsigned int number=-1
        ) -> int {
            if (mode == NUMBER) {
                return 2;
            }
            if (!merged[vertice].isHiddenInfo()) {
                return merged[vertice].leaf;
            }
            return trees[number][merged[vertice].getHiddenInfo(number)].leaf;
        }
    );
    vector <unsigned int> output;
    vector <IndexedPair<unsigned int> > lcaSuffixes;
    findLcaSuffix(merged.root, merged, trees, lcaSuffixes);

    vector <unsigned int> toVisit;
    toVisit.push_back(merged.root);

    while (toVisit.size()) {
        unsigned int v = toVisit.back();
        toVisit.pop_back();
        if (merged[v].isHiddenInfo()
            && (output.size() <= v || output[v] == -1)) {
            trueLengthFillDfs(v, merged, getter, lcaSuffixes[v],
                inputLength, output
            );
        }
        for (auto const &u: merged[v]) {
            toVisit.push_back(u);
        }
    }
    return output;
}

void correctMerge(unsigned int v, unsigned int parentsPlace, SuffixTree &merged,
    IndexedPair<MergeTreesStruct> &updatedTrees,
    const vector <unsigned int> &trueLength, const vector <int> &input,
    unsigned int copyTree
) {
    if (merged[v].isHiddenInfo()) {
        doSomething(updatedTrees, [&merged, &v] (MergeTreesStruct &tree) {
            tree.evaluate(merged, v);
        });
        if (copyTree == 2 && merged[v].depth != trueLength[v]) {
            unsigned int commonLength = trueLength[v]
                - merged[merged[v].parent].depth;
            MergeTreesStruct& donor = minimal(updatedTrees,
                [&input, &commonLength] (MergeTreesStruct &tree)
                    -> int {
                    auto const &node = tree.tree[tree.info];
                    return input[tree.suffix[tree.info]
                        + (node.parent == -1 ? 0 : tree.tree[node.parent].depth)
                        + commonLength];
                }
            );
            copyTree = donor.number;
            copyNodeExceptParentAndChildren(donor.tree[donor.info], merged[v]);
            unsigned int newNodeIndex = merged.splitEdge(merged[v].parent,
                parentsPlace, commonLength
            );
            unsigned int newCopy = merged.newNode(newNodeIndex);
            merged[newNodeIndex].push_back(newCopy);
            MergeTreesStruct &notDonor = *xorPointers(&donor, &updatedTrees[0],
                &updatedTrees[1]
            );
            copyNodeExceptParentAndChildren(notDonor.tree[notDonor.info],
                merged[newCopy]
            );
            merged[newCopy].indexOfParentEdge += commonLength;
            copySubTree(notDonor.tree, merged, notDonor.info, newCopy);
        } else {
            MergeTreesStruct &donor = minimal(updatedTrees,
                [&copyTree] (MergeTreesStruct &tree) -> pair<bool, bool> {
                    return make_pair(copyTree == 2 || copyTree != tree.number,
                        !(tree.tree[tree.info].leaf != -1)
                    );
                }
            );
            copyNodeExceptParentAndChildren(donor.tree[donor.info], merged[v]);
        }
    }
    if (copyTree != 2 && merged[v].leaf != -1
            && merged[v].leaf % 2 != copyTree) {
        merged[v].leaf = -1;
    }
    for (unsigned int i = 0; i < merged[v].size(); ++i) {
        unsigned int u = merged[v][i];
        correctMerge(u, i, merged, updatedTrees, trueLength, input, copyTree);
    }
}

void checkTree(const char *message, const SuffixTree &tree,
    const vector <int> &input
) {
    #ifdef _DEBUG
        #ifdef _PRINT_DBG
            printf(message);
        #endif
        checkDfs(tree, tree.root, input);
    #endif
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

    checkTree("NEW TREE1:\n", tree1, input);
    checkTree("NEW TREE2:\n", tree2, input);

    IndexedPair<MergeTreesStruct> updatedTrees(MergeTreesStruct(tree1, 0),
        MergeTreesStruct(tree2, 1)
    );

    vector <unsigned int> trueLength = computeTrueLength(merged,
        IndexedPair<const SuffixTree &>(tree1, tree2), input.size()
    );

    correctMerge(merged.root, -1, merged, updatedTrees, trueLength, input);
    return merged;
}

int cleanTreeDfs(unsigned int v, unsigned int parent, SuffixTree &tree) {
    int leaf = tree[v].leaf;
    for (unsigned int i = 0; i < tree[v].size(); ++i) {
        unsigned int u = tree[v][i];
        int newLeaf = cleanTreeDfs(u, i, tree);
        if (leaf == -1) {
            leaf = newLeaf;
        }
    }
    unsigned int freeCell = 0;
    for (unsigned int i = 0; i < tree[v].size(); ++i) {
        if (tree[v][i] != -1) {
            tree[v][freeCell] = tree[v][i];
            ++freeCell;
        }
    }
    tree[v].resize(freeCell);
    if (tree[v].size() <= 1 && tree[v].leaf == -1 && v != tree.root) {
        tree.deleteUselessNode(v, parent, leaf);
    }
    return leaf;
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
    checkTree("EVEN:\n", even, input);

    SuffixTree odd = buildOddSuffixTree(even, input);
    checkTree("ODD:\n", odd, input);

    SuffixTree almostResult = mergeTrees(even, odd, input);
    checkTree("ALMOST:\n", almostResult, input);

    cleanTreeDfs(almostResult.root, -1, almostResult);
    checkTree("RESL:\n", almostResult, input);

    return almostResult;
}

SuffixTree buildSuffixTree(const vector <int> &input) {
    #ifdef _PRINT_DBG
        for (unsigned int i = 0; i < input.size(); ++i) {
            printf("%d ", input[i]);
        }
        printf("\n");
    #endif
    SuffixTree result = buildTempSuffixTree(input);
    result[result.root].leaf = input.size();
    return result;
}

unsigned long long countSubstrings(SuffixTree &tree, unsigned int v) {
    tree.checkNode(tree[v]);
    unsigned long long sum = tree[v].lengthOfEdge(tree);
    for (auto const &u: tree[v]) {
        sum += countSubstrings(tree, u);
    }
    return sum;
}

void randGen(int length, int alph, int tests) {
    static vector <int> res;
    for (int i = 0; i < tests; ++i) {
        res.resize(length);
        for (int j = 0; j < length; ++j) {
            res[j] = rand() % alph;
            printf("%d ", res[j]);
        }
        printf("\n");
        buildSuffixTree(res);
    }
}

int main() {
    // randGen(10, 3, 100);
    std::string s;
    std::cin >> s;

    unsigned int n = s.size();
    vector <int> input(n);
    for (unsigned int i = 0; i < n; ++i) {
        input[i] = s[i] - 'a';
    }
    SuffixTree tree = buildSuffixTree(input);
    std::cout << countSubstrings(tree, tree.root) << std::endl;
    return 0;
}
