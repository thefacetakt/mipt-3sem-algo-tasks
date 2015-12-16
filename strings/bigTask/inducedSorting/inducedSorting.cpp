#include <cstdio>
#include <vector>
#include <algorithm>
#include <climits>
#include <ctime>
#include "inducedSorting.hpp"


using std::vector;
using std::min;
using std::max;


vector <EType> detectTypes(const vector <unsigned int> &input,
    vector <unsigned int> &lmsPositions) {

    vector <EType> type(input.size());

    type.back() = S;
    for (int i = input.size() - 2; i >= 0; --i) {
        if (input[i] < input[i + 1]) {
            type[i] = S;
        } else if (input[i] > input[i + 1]) {
            type[i] = L;
        } else {
            type[i] = type[i + 1];
        }
        if (type[i] == L && type[i + 1] == S) {
            type[i + 1] = LMS;
            lmsPositions.push_back(i + 1);
        }
    }
    type.back() = LMS;
    reverse(lmsPositions.begin(), lmsPositions.end());
    return type;
}

vector <unsigned int> generateReducedAlphabet(
    const vector <unsigned int> &input,
    const vector <EType> &type,
    const vector <unsigned int> &suffixArray) {

    vector<unsigned int> newAlphabet(input.size(), -1);

    unsigned int lastLms = -1;
    unsigned int lmsCount = 0;
    for (const auto &index: suffixArray) {
        if (type[index] == LMS) {
            ++lmsCount;
            if (lastLms == -1) {
                newAlphabet[index] = 0;
            } else {
                unsigned int i;
                for (i = 1; type[lastLms + i] != LMS
                    && input[lastLms + i] == input[index + i]; ++i) {
                }
                newAlphabet[index] = newAlphabet[lastLms];
                if (input[lastLms] != input[index]
                    || (type[lastLms + i] != LMS || type[index + i] != LMS)) {
                    ++newAlphabet[index];
                }
            }
            lastLms = index;
        }
    }
    return newAlphabet;
}

vector <unsigned int> generateReducedString(
    const vector <unsigned int> &newAlphabet,
    const vector <unsigned int> &lmsPositions) {

    vector <unsigned int> newInput;
    newInput.reserve(lmsPositions.size());
    for (unsigned int index: lmsPositions) {
        newInput.push_back(newAlphabet[index]);
    }
    return newInput;
}

void induce(const vector <unsigned int> &input,
    const vector <EType> &type,
    vector <unsigned int> &suffixArray,
    vector <unsigned int> tail,
    unsigned int alphabetSize) {

    vector <unsigned int> head(alphabetSize, 0);

    for (unsigned int i = 1; i < tail.size(); ++i) {
        head[i] = tail[i - 1];
    }

    for (auto const &index: suffixArray) {
        if (index != -1 && index != 0 &&
            type[index - 1] == L) {
            suffixArray[head[input[index - 1]]++] = index - 1;
        }
    }


    for (int i = suffixArray.size() - 1; i >= 0; --i) {
        unsigned int index = suffixArray[i];
        if (index != -1 && index != 0 &&
            (type[index - 1] == S || type[index - 1] == LMS)) {
            suffixArray[--tail[input[index - 1]]] = index - 1;
        }
    }
}

vector <unsigned int> countTail(const vector <unsigned int> &input,
    unsigned int alphabetSize) {

    vector <unsigned int> tail(alphabetSize, 0);
    for (const auto &c: input) {
        ++tail[c];
    }
    for (unsigned int i = 1; i < alphabetSize; ++i) {
        tail[i] += tail[i - 1];
    }

    return tail;
}


vector <unsigned int> inducedSortingPrepared(const vector <unsigned int> &input,
    int alphabetSize) {

    if (input.size() == 1) {
        return vector <unsigned int> (1, 0);
    }

    vector <unsigned int> lmsPositions;
    vector <EType> type = detectTypes(input, lmsPositions);


    vector <unsigned int> tail = countTail(input, alphabetSize);
    vector <unsigned int> tailBackUp = tail;

    vector <unsigned int> suffixArray(input.size(), -1);

    for (const auto &index: lmsPositions) {
        suffixArray[--tail[input[index]]] = index;
    }

    induce(input, type, suffixArray, tailBackUp, alphabetSize);

    vector <unsigned int> reducedAlphabet
        = generateReducedAlphabet(input, type, suffixArray);

    vector <unsigned int> reducedString
        = generateReducedString(reducedAlphabet, lmsPositions);

    vector <unsigned int> reducedSuffixArray
        = inducedSortingPrepared(reducedString, reducedAlphabet.size());

    tail = tailBackUp;
    suffixArray.assign(input.size(), -1);

    for (int i = reducedSuffixArray.size() - 1; i >= 0; --i) {
        unsigned int index = reducedSuffixArray[i];
        suffixArray[--tail[input[lmsPositions[index]]]] = lmsPositions[index];
    }

    induce(input, type, suffixArray, tailBackUp, alphabetSize);

    return suffixArray;
}

vector <unsigned int> inducedSorting(vector <unsigned int> input) {
    unsigned int alphabetSize = 0;
    for (auto &c: input) {
        ++c;
        alphabetSize = max(c, alphabetSize);
    }
    ++alphabetSize;
    input.push_back(0);

    vector <unsigned int> answer = inducedSortingPrepared(input, alphabetSize);

    return vector <unsigned int> (answer.begin() + 1, answer.end());
}


vector <unsigned int> gen(int n) {
    vector <unsigned int> res;
    for (int i = 0; i < n; ++i) {
        res.push_back(rand() % 26 + 'a');
    }
    return res;
}

int main() {
    vector <unsigned int> s = gen(1e7);
    vector <unsigned int> suffixArray = inducedSorting(s);

    printf("%.6lf\n", (clock() * 1.0 / CLOCKS_PER_SEC));
    unsigned long long ans = 0;
    for (const auto &i: suffixArray) {
        ans += i;
    }
    printf("%llu\n", ans);

    return 0;
}
