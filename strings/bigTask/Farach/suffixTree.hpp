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
