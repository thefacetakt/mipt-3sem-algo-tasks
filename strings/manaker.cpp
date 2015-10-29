#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>
#include <utility>

const size_t MAXN = 500000;
char s[MAXN];
int d[MAXN];

int main() {
    
    scanf("%s", s);
    int n = strlen(s);
    
    
    int answerLeft = 0;
    int answerRight = 0;
    
    for (int z = 0; z < 2; ++z) {
        int left = 0, right = -1;
        for (int i = 0; i < n; ++i) {
            if (i > right) {
                d[i] = 1 - z;
            } else {
                d[i] = std::min(d[left + right - i + z], right - i + z);
            }
            while (d[i] + z <= i && i + d[i] < n && s[i + d[i]] == s[i - d[i] - z]) {
                ++d[i];
            }
            if (i + d[i] - z > right) {
                left = i - d[i] + (1 - z);
                right = i + d[i] - 1;
            }
            if (answerRight - answerLeft < i + d[i] - (i - d[i] - z) - 1) {
                answerLeft = i - d[i] - z + 1;
                answerRight = i + d[i];
            }
        }
    }
    
    s[answerRight] = '\0';
    printf("%s\n", s + answerLeft);
    return 0;
         
}