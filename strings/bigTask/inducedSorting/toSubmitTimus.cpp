#include <climits>

#ifndef _USEFUL_STRUCTURES_IS
#define _USEFUL_STRUCTURES_IS

#include <vector>

using std::vector;

enum Type {
    MINUS,
    PLUS,
    STAR,
    TOTAL
};

struct Node {
    vector <unsigned int> children;
    vector <unsigned int> term;
};

struct StringSortingItem {
    unsigned int position;
    unsigned int string;

    StringSortingItem();

    StringSortingItem(unsigned int position, unsigned int string);
};

#endif



StringSortingItem::StringSortingItem() {
}

StringSortingItem::StringSortingItem(unsigned int position, unsigned int string)
    : position(position), string(string) {
}

#ifndef _INDUCED_SORTING
#define _INDUCED_SORTING

#include <vector>
#include <algorithm>


using std::vector;
using std::max;

template<class T, class Compare>
vector <T> radixSort(const vector <T> &input, Compare comp) {
    unsigned int maxElement = 0;
    for (auto const &element: input) {
        maxElement = max(maxElement, static_cast<unsigned int>(comp(element)));
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

vector <unsigned int> finalDfs(const vector <vector <unsigned int> > &strings,
    const vector<Node> &trie);

void addNode(const StringSortingItem &item, vector <Node> &trie,
    vector<unsigned int> &position,
    const vector <vector <unsigned int> > &strings);

vector <unsigned int> sortOfString(
    const vector <vector <unsigned int> > &strings);

vector <Type> detectTypes(const vector <unsigned int> &input);

void sortStars(const vector <unsigned int> &input,
    unsigned int maxSymbol,
    const vector <Type> &type,
    vector <vector <unsigned int> > &out);

void induceMinuses(const vector <unsigned int> &input,
    const vector <Type> &type,
    vector <vector <vector <unsigned int> > > &sortedParts);

void inducePluses(const vector <unsigned int> &input, const vector<Type> &type,
    vector <vector <vector <unsigned int> > > &sortedParts);

unsigned int defineMaxSymbolAndAddZeroSymbol(vector <unsigned int> &input);

vector <unsigned int> restoreAnswer(unsigned int maxSymbol,
    vector <vector <vector <unsigned int> > > &sortedParts);


vector <unsigned int> inducedSortingChangable(vector <unsigned int> &input);

vector <unsigned int> inducedSorting (vector <unsigned int> input);

#endif

#include <cstdio>
#include <iostream>
#include <vector>
#include <algorithm>
#include <climits>
#include <utility>

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

#include <algorithm>

using std::max;
using std::min;

vector <unsigned int> lcpArray(const vector <unsigned int> &input,
    const vector <unsigned int> &suffixArray
) {
    vector <unsigned int> antiArray(suffixArray.size());
    for (unsigned int i = 0; i < suffixArray.size(); ++i) {
        antiArray[suffixArray[i]] = i;
    }
    vector <unsigned int> lcp(suffixArray.size() - 1, 0);
    unsigned int current = 0;
    for (unsigned int i = 0; i < suffixArray.size(); ++i) {
        if(antiArray[i] == suffixArray.size() - 1) {
                current = 0;
                continue;
        }
        lcp[antiArray[i]] = current;
        for (unsigned int i1 = i + current,
            i2 = suffixArray[antiArray[i] + 1] + current;
            i1 < input.size() && i2 < input.size() && input[i1] == input[i2];
            ++i1, ++i2, ++lcp[antiArray[i]]) {
        }
        current = max(static_cast<int>(lcp[antiArray[i]]) - 1, 0);
    }
    return lcp;
}

unsigned long long countSubstrings(const vector <unsigned int> &lcp,
    const vector <unsigned int> &suffixArray,
    unsigned int length
) {
    unsigned long long answer = length  - suffixArray[0];
    for (unsigned int i = 1; i < length; ++i) {
        answer += (length - suffixArray[i] - lcp[i - 1]);
    }
    return answer;
}

unsigned long long countSubstringsComplete(const vector <unsigned int> &input) {
    vector <unsigned int> suffixArray = inducedSorting(input);
    return countSubstrings(lcpArray(input, suffixArray),
        suffixArray, input.size()
    );
}

void relax(unsigned int &maximum, unsigned int k, unsigned int j,
    unsigned int i, unsigned int &currentLcp,
    const vector <unsigned int> &suffix,
    const vector <unsigned int> &lcp
) {
    if (i < suffix[j] && suffix[j] < i + k) {
        maximum = max(maximum, min(i + k - suffix[j], currentLcp));
    }
}

unsigned int maxStart(unsigned int i, unsigned int k,
    const vector <unsigned int> &lcp, const vector <unsigned int> &suffix
) {
    unsigned int place = 0;
    for (place = 0; place < suffix.size(); ++place) {
        if (suffix[place] == i) {
            break;
        }
    }
    unsigned int maximum = 0;
    unsigned int currentLcp = lcp[place];
    for (unsigned int j = place + 1; j < suffix.size(); ++j) {
        relax(maximum, k, j, i, currentLcp, suffix, lcp);
        if (j + 1 != suffix.size()) {
            currentLcp = min(currentLcp, lcp[j]);
        }
    }

    currentLcp = suffix.size() - i;
    for (unsigned int j = place - 1; j != UINT_MAX; --j) {
        currentLcp = min(currentLcp, lcp[j]);
        relax(maximum, k, j, i, currentLcp, suffix, lcp);
    }
    return k - min(k, maximum);;
}

int main() {
    std::string s;
    unsigned int k;
    std::cin >> k >> s;

    unsigned int n = s.size();
    s += s;
    vector <unsigned int> input(k);

    for (unsigned int i = 0; i < k; ++i) {
        input[i] = s[i] - 'a';
    }
    unsigned long long start = countSubstringsComplete(input);
    printf("%llu ", start);
    input.resize(2 * n);

    for (unsigned int i = k; i < 2 * n; ++i) {
        input[i] = s[i % n] - 'a';
    }

    vector <unsigned int> reverseInput(input.rbegin(), input.rend());
    vector <unsigned int> suffixArray = inducedSorting(input),
        reverseArray = inducedSorting(reverseInput);
    vector <unsigned int> lcp = lcpArray(input, suffixArray),
        reverseLcp = lcpArray(reverseInput, reverseArray);

    for (unsigned int i = 1; i < n; ++i) {
        start -= maxStart(i - 1, k, lcp, suffixArray);
        start += maxStart(2 * n - i - k, k, reverseLcp, reverseArray);
        printf("%llu ", start);
    }
    printf("\n");

    return 0;
}
