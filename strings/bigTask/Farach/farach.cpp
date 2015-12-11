#include <cstdio>
#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>
#include <climits>
#include <cassert>
#include <functional>
#include <string>
#include "eulerPair.hpp"
#include "LCA.hpp"
#include "suffixTree.hpp"
#include "usefulStructures.hpp"
#include "farach.hpp"


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
