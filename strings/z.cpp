#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

struct TwoStrings {
    std::string s, t;
    static char sentinel;
    
    TwoStrings() {
    }
    
    char& operator[](size_t i) {
        if (i < s.size()) {
            return s[i];
        }
        if (i == s.size()) {
            return sentinel;
        };
        return t[i - s.size() - 1];
    }
    
    size_t size() {
        return s.size() + t.size() + 1;
    }
};
char TwoStrings::sentinel = '#';

std::istream& operator>>(std::istream &in, TwoStrings &strings) {
    in >> strings.s >> strings.t;
    return in;
}

int main() {
    TwoStrings st;
    std::cin >> st;
    
    std::vector<size_t> z(st.s.size(), 0);
    size_t left = 0, right = 0;
    
    for (size_t i = 1; i < st.size(); ++i) {
        size_t currentZ = 0;
        if (i <= right) {
            currentZ = std::min(right - i + 1, z[i - left]);
        }
        while (i + currentZ < st.size() && st[currentZ] == st[i + currentZ]) {
            ++currentZ;
        }
        if (i + currentZ - 1 > right) {
            left = i;
            right = i + currentZ - 1;
        }
        if (i < st.s.size()) {
            z[i] = currentZ;
        }
        if (i >st.s.size()) {
            if (currentZ == st.s.size()) {
                std::cout << i - st.s.size() - 1 << " ";
            }
        }
    }
    std::cout << "\n";
    return 0;
}