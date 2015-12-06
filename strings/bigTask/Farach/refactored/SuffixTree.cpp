#include <vector>
#include <cassert>
#include "SuffixTree.hpp"

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
