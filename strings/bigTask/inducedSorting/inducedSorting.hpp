#include <cstdio>
#include <vector>
#include <algorithm>
#include <climits>
#include <ctime>

using std::vector;
using std::min;
using std::max;

enum EType {
    L,
    S,
    LMS
};

vector <EType> detectTypes(const vector <unsigned int> &input,
    vector <unsigned int> &lmsPositions);

vector <unsigned int> generateReducedAlphabet(
    const vector <unsigned int> &input,
    const vector <EType> &type,
    const vector <unsigned int> &suffixArray);

vector <unsigned int> generateReducedString(
    const vector <unsigned int> &newAlphabet,
    const vector <unsigned int> &lmsPositions);

void induce(const vector <unsigned int> &input,
    const vector <EType> &type,
    vector <unsigned int> &suffixArray,
    vector <unsigned int> tail,
    unsigned int alphabetSize);

vector <unsigned int> countTail(const vector <unsigned int> &input,
    unsigned int alphabetSize);


vector <unsigned int> inducedSortingPrepared(const vector <unsigned int> &input,
    int alphabetSize);

vector <unsigned int> inducedSorting(vector <unsigned int> input);

vector <unsigned int> lcpArray(const vector <unsigned int> &input,
    const vector <unsigned int> &suffixArray);

unsigned long long countSubstrings(const vector <unsigned int> &lcp,
    const vector <unsigned int> &suffixArray,
    unsigned int length);

vector <unsigned int> gen(int n);

int main();
