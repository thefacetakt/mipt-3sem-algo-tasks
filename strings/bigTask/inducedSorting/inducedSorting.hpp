#ifndef _INDUCED_SORTING
#define _INDUCED_SORTING

#include <vector>
#include <algorithm>
#include "usefulStructures.hpp"

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
