#include <iostream>
#include <string>
#include <vector>

int main() {
    std::string s[2];
    size_t n[2];
    std::vector <size_t> pi;
    
    for (size_t z = 0; z < 2; ++z) {
        std::cin >> s[z];
        n[z] = s[z].size();
        if(z == 0) {
            pi.resize(n[z]);
            pi[0] = 0;
        }
        size_t currentPi = 0;
        for (size_t i = 1 - z; i < n[z]; ++i) {
            size_t newPi = currentPi;
            while (newPi > 0 && s[0][newPi] != s[z][i]) {
                newPi = pi[newPi - 1];
            }
            if (s[0][newPi] == s[z][i]) {
                ++newPi;
            }
            if (z == 0) {
                pi[i] = newPi;
            } else {
                if (newPi == n[0]) {
                    std::cout << i - n[0] + 1 << " ";
                }
            }
            currentPi = newPi;
        }
    }
    std::cout << "\n";
}