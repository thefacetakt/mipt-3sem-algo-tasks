#include <cstdio>
#include <vector>
#include <algorithm>
#include <climits>

using std::vector;
using std::min;
using std::max;

enum EType {
    L,
    S,
    LMS
};

vector <EType> detectTypes(const vector <unsigned int> &input,
    vector <unsigned int> &lmsPositions) {

    vector <EType> type(input.size());

    type.back() = S;
    for (unsigned int i = input.size() - 2; i != UINT_MAX; --i) {
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

template<typename T>
void print(const vector <T> &x) {
    for (auto const &i: x) {
        printf("%d ", i);
    }
    printf("\n");
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


    for (unsigned int i = suffixArray.size() -1; i != UINT_MAX; --i) {
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
    // for (int i = 0; i < input.size(); ++i) {
    //     if (input[i] != 0) {
    //         printf("%c", input[i] + 'a' - 1);
    //     }
    // }
    // printf("\n");
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

    for (unsigned int i = reducedSuffixArray.size() - 1; i != UINT_MAX; --i) {
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

#include <string>
#include <iostream>
using std::string;
using std::cin;
using std::cout;

vector <unsigned int> lcpArray(const vector <unsigned int> &input,
    const vector <unsigned int> &suffixArray) {
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
    unsigned int length) {
    unsigned long long answer = length  - suffixArray[0];
    for (unsigned int i = 1; i < length; ++i) {
        answer += (length - suffixArray[i] - lcp[i - 1]);
    }
    return answer;
}

std::string gen(int n) {
    std::string res;
    for (int i = 0; i < n; ++i) {
        res += (rand() % 26 + 'a');
    }
    return res;
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
    const vector <unsigned int> &lcp) {
    if (i < suffix[j] && suffix[j] < i + k) {
        maximum = max(maximum, min(i + k - suffix[j], currentLcp));
    }
}

unsigned int maxStart(unsigned int i, unsigned int k,
    const vector <unsigned int> &lcp, const vector <unsigned int> &suffix) {
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
    vector <unsigned int> input(2 * n);

    for (unsigned int i = 0; i < 2 * n; ++i) {
        input[i] = s[i] - 'a';
    }
    for (int i = 0; i < n; ++i) {
        printf("%llu ", countSubstringsComplete(
            vector<unsigned int> (input.begin() + i, input.begin() + i + k)));
    }
    printf("\n");

    return 0;
}
