#include <cstdio>
#include <iostream>
#include <vector>
#include <algorithm>
#include <climits>
#include <utility>
#include "inducedSorting.hpp"
#include "usefulStructures.hpp"

using std::vector;
using std::pair;
using std::max;

vector <unsigned int> finalDfs(const vector <vector <unsigned int> > &strings,
    const vector<Node> &trie
) {
    vector <unsigned int> ans(strings.size());
    vector <unsigned int> tempAns;
    vector <unsigned int> stack;
    stack.push_back(0);
    unsigned int currentNumber = 0;
    while (stack.size()) {
        unsigned int v = stack.back();
        stack.pop_back();
        for (auto const &i: trie[v].term) {
            ans[i] = currentNumber;
            tempAns.push_back(i);
        }
        if (trie[v].term.size()) {
            ++currentNumber;
        }
        for (unsigned int i = trie[v].children.size() - 1; i != UINT_MAX; --i) {
            stack.push_back(trie[v].children[i]);
        }

    }
    return ans;
}

void addNode(const StringSortingItem &item, vector <Node> &trie,
    vector<unsigned int> &position,
    const vector <vector <unsigned int> > &strings
) {
    Node &node = trie[position[item.string]];
    if (!node.children.size()
            || strings[trie[node.children.back()].term[0]][item.position]
            != strings[item.string][item.position]
        ) {
        unsigned int newNodeIndex = trie.size();
        trie.push_back(Node());
        trie[position[item.string]].children.push_back(newNodeIndex);
        trie.back().term.push_back(item.string);
        position[item.string] = newNodeIndex;
    } else {
        position[item.string] = node.children.back();
        trie[node.children.back()].term.push_back(item.string);
    }
}

vector <unsigned int> sortOfString(
    const vector <vector <unsigned int> > &strings
) {
    vector <StringSortingItem> temp;
    for (unsigned int i = 0; i < strings.size(); ++i) {
        for (unsigned int j = 0; j < strings[i].size(); ++j) {
            temp.push_back(StringSortingItem(j, i));
        }
    }

    temp = radixSort(radixSort(temp,
        [&strings] (const StringSortingItem &item) -> unsigned int {
            return strings[item.string][item.position];
        }), [] (const StringSortingItem &item) -> unsigned int {
            return item.position;
        }
    );

    vector <Node> trie(1);
    vector <unsigned int> position(strings.size());
    for (unsigned int i = 0; i < strings.size(); ++i) {
        position[i] = 0;
        trie[0].term.push_back(i);
    }

    for (auto const &item: temp) {
        addNode(item, trie, position, strings);
    }

    for (auto &node: trie) {
        node.term.clear();
    }
    for (unsigned int i = 0; i < position.size(); ++i) {
        trie[position[i]].term.push_back(i);
    }

    return finalDfs(strings, trie);
}

vector <Type> detectTypes(const vector <unsigned int> &input) {
    vector <Type> type(input.size());
    for (unsigned int i = input.size() - 1; i != UINT_MAX; --i) {
        if (i + 1 == input.size()) {
            type[i] = MINUS;
        } else if (type[i + 1] == MINUS) {
            if (input[i] >= input[i + 1]) {
                type[i] = MINUS;
            } else {
                type[i] = PLUS;
            }
        } else {
            if (input[i] <= input[i + 1]) {
                type[i] = PLUS;
            } else {
                type[i] = MINUS;
                type[i + 1] = STAR;
            }
        }
    }
    return type;
}

void sortStars(const vector <unsigned int> &input,
    unsigned int maxSymbol,
    const vector <Type> &type,
    vector <vector <unsigned int> > &out
) {
    vector <vector <unsigned int> > starStrings;
    vector <unsigned int> currentString;
    for (unsigned int i = 0; i < input.size(); ++i) {
        if (i + 1 == input.size() || type[i] == STAR) {
            if (currentString.size()) {
                currentString.push_back(input[i]);
                currentString.push_back(maxSymbol);
                starStrings.push_back(currentString);
                currentString.clear();
            }
            currentString.push_back(input[i]);
        } else {
            if (currentString.size()) {
                currentString.push_back(input[i]);
            }
        }
    }
    vector <unsigned int> newSymbols = sortOfString(starStrings);
    vector <unsigned int> newString;

    vector <unsigned int> starPositions;
    for (unsigned int i = 0; i < input.size(); ++i) {
        if (type[i] == STAR) {
            newString.push_back(newSymbols[starPositions.size()]);
            starPositions.push_back(i);
        }
    }

    vector <unsigned int> sortedStars = inducedSortingChangable(newString);
    for (unsigned int i = 0; i < sortedStars.size(); ++i) {
        out[input[starPositions[sortedStars[i]]]].push_back(
            starPositions[sortedStars[i]]
        );
    }
}

void induceMinuses(const vector <unsigned int> &input,
    const vector <Type> &type,
    vector <vector <vector <unsigned int> > > &sortedParts
) {
    sortedParts[MINUS][0].push_back(input.size() - 1);
    vector <unsigned int> ptr(sortedParts[MINUS].size(), 0);
    for (unsigned int i = 0; i < ptr.size(); ++i) {
        while (ptr[i] != sortedParts[MINUS][i].size()) {
            unsigned int index = sortedParts[MINUS][i][ptr[i]];
            if (index > 0 && type[index - 1] == MINUS) {
                sortedParts[MINUS][input[index - 1]].push_back(index - 1);
            }
            ++ptr[i];
        }
        for (unsigned int j = 0; j  < sortedParts[STAR][i].size(); ++j) {
            sortedParts[MINUS][input[sortedParts[STAR][i][j] - 1]].push_back(sortedParts[STAR][i][j] - 1);
        }
    }
}

void inducePluses(const vector <unsigned int> &input, const vector<Type> &type,
    vector <vector <vector <unsigned int> > > &sortedParts
) {
    vector <unsigned int> ptr(sortedParts[PLUS].size(), 0);
    for (unsigned int i = ptr.size() - 1; i != UINT_MAX; --i) {
        while (ptr[i] != sortedParts[PLUS][i].size()) {
            unsigned int index = sortedParts[PLUS][i][ptr[i]];
            if (index > 0 && type[index - 1] != MINUS) {
                sortedParts[PLUS][input[index - 1]].push_back(index - 1);
            }
            ++ptr[i];
        }
        for (unsigned int j = sortedParts[MINUS][i].size() - 1;
            j != UINT_MAX; --j
        ) {
            int index = sortedParts[MINUS][i][j];
            if (index > 0 && type[index - 1] != MINUS) {
                sortedParts[PLUS][input[index - 1]].push_back(index - 1);
            }
        }
    }
}

unsigned int defineMaxSymbolAndAddZeroSymbol(vector <unsigned int> &input) {
    unsigned int maxSymbol = 0;
    for (int i = 0; i < input.size(); ++i) {
        input[i] += 1;
        maxSymbol = max(maxSymbol, input[i] + 1);
    }
    input.push_back(0);
    return maxSymbol;
}

vector <unsigned int> restoreAnswer(unsigned int maxSymbol,
    vector <vector <vector <unsigned int> > > &sortedParts
) {
    vector <unsigned int> answer;
    for (unsigned int i = 1; i < maxSymbol; ++i) {
        for (auto const &j: sortedParts[MINUS][i]) {
            answer.push_back(j);
        }
        for (unsigned int j = sortedParts[PLUS][i].size() - 1;
            j != UINT_MAX; --j
        ) {
            answer.push_back(sortedParts[PLUS][i][j]);
        }
    }
    return answer;
}

vector <unsigned int> inducedSortingChangable(vector <unsigned int> &input) {
    if (input.size() == 0) {
        return vector <unsigned int> ();
    }
    if (input.size() == 1) {
        vector <unsigned int> ans(1, 0);
        return ans;
    }

    unsigned int maxSymbol = defineMaxSymbolAndAddZeroSymbol(input);

    vector <Type> type = detectTypes(input);
    vector <vector <vector <unsigned int> > > sortedParts(
        static_cast<unsigned int>(TOTAL),
        vector<vector<unsigned int> > (maxSymbol)
    );

    sortStars(input, maxSymbol, type, sortedParts[STAR]);
    induceMinuses(input, type, sortedParts);
    inducePluses(input, type, sortedParts);


    return restoreAnswer(maxSymbol, sortedParts);
}

vector <unsigned int> inducedSorting (vector <unsigned int> input) {
    return inducedSortingChangable(input);
}
