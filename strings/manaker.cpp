#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>
#include <utility>

const size_t MAXN = 500000;
char s[MAXN];
int radius[MAXN];

int main() {
    
    scanf("%s", s);
    int n = strlen(s);
    
    
    int answerLeft = 0;
    int answerRight = 0;
    
    for (int z = 0; z < 2; ++z) {
        int left = 0, right = -1;
        for (int i = 0; i < n; ++i) {
            if (i > right) {
                radius[i] = 1 - z;
            } else {
                radius[i] = std::min(radius[left + right - i + z], right - i + z);
            }
            while (radius[i] + z <= i && i + radius[i] < n && s[i + radius[i]] == s[i - radius[i] - z]) {
                ++radius[i];
            }
            if (i + radius[i] - z > right) {
                left = i - radius[i] + (1 - z);
                right = i + radius[i] - 1;
            }
            if (answerRight - answerLeft < i + radius[i] - (i - radius[i] - z) - 1) {
                answerLeft = i - radius[i] - z + 1;
                answerRight = i + radius[i];
            }
        }
    }
    
    s[answerRight] = '\0';
    printf("%s\n", s + answerLeft);
    return 0;
         
}