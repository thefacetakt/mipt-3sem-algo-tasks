#define DEBUG

#include <cstdio>
#include <vector>
#include <algorithm>
#include <cstdarg>
#include <cassert>

using std::vector;
using std::max;
using std::min;
using std::swap;
using std::abs;

const unsigned int MAXN = 128;
const unsigned int MAXN_1 = 127;
const unsigned int LOG_MAXN = 7;


const int WHITE = 0;
const int BLACK = 1;
const int EMPTY = 2;
const int OUTER = 4;


const unsigned int MAX_BUFFER = 100;

const int W = 8, H = 8;

int map[W + 2][H + 2];
int availible[W + 2][H + 2];

int directionsX[W];
int directionsY[H];


struct FastPositionsStack {
    int antiPos[W + 5][H + 5];
    int xPos[MAXN];
    int yPos[MAXN];
    int ptr;

    FastPositionsStack() {
        ptr = 0;
    }

    void push(int x, int y) {
        antiPos[x][y] = ptr;
        xPos[ptr] = x;
        yPos[ptr] = y;
        ++ptr;
    }

    void pop() {
        --ptr;
    }

    int topX() {
        #ifdef DEBUG
            assert(!empty());
        #endif
        return xPos[ptr - 1];
    }

    int topY() {
        #ifdef DEBUG
            assert(!empty());
        #endif
        return yPos[ptr - 1];
    }

    int top(int what) {
        #ifdef DEBUG
            assert(0 <= what && what <= 1);
        #endif
        return (what == 0 ? topX() : topY());
    }

    void popElement(int pos) {
        #ifdef DEBUG
            assert(0 <= pos && pos < ptr);
        #endif
        antiPos[xPos[ptr - 1]][yPos[ptr - 1]] = pos;
        swap(xPos[pos], xPos[ptr - 1]);
        swap(yPos[pos], yPos[ptr - 1]);
        pop();
    }

    bool empty() {
        return ptr == 0;
    }
};

struct FastPositionsQueue {
    int xPos[MAXN];
    int yPos[MAXN];
    int begin;
    int end;

    FastPositionsQueue() {
        begin = end = 0;
    }

    void push(int x, int y) {
        xPos[end] = x;
        yPos[end] = y;
        end = (end + 1) & MAXN_1;
    }

    void pop() {
        begin = (begin + 1) & MAXN_1;
    }

    int topX() {
        #ifdef DEBUG
            assert(!empty());
        #endif
        return xPos[begin];
    }

    int topY() {
        #ifdef DEBUG
            assert(!empty());
        #endif
        return yPos[begin];
    }

    int top(int what) {
        #ifdef DEBUG
            assert(0 <= what && what <= 1);
        #endif
        return (what == 0 ? topX() : topY());
    }

    bool empty() {
        return begin == end;
    }
};

bool myColor;

FastPositionsStack myStack;


void fprint(std::FILE * stream, const char* format, ...) {
    va_list args;
    va_start(args, format);
    vfprintf(stream, format, args);
    fflush(stream);
    va_end(args);
}

void eprint(const char* format, ...) {
    va_list args;
    va_start(args, format);
    vfprintf(stderr, format, args);
    fflush(stderr);
    va_end(args);
}

void print(const char* format, ...) {
    va_list args;
    va_start(args, format);
    vprintf(format, args);
    fflush(stdout);
    va_end(args);
}

char charState(int state) {
    if (state == WHITE) {
        return 'w';
    }
    if (state == BLACK) {
        return 'b';
    }
    if (state == EMPTY) {
        return '.';
    }
    return '#';
}


inline void printBoard() {
    #ifdef DEBUG
        for (int i = 0; i < W + 2; ++i) {
            for (int j = 0; j < H + 2; ++j) {
                eprint("%c", charState(map[i][j]));
            }
            eprint("\n");
        }
    #endif
}

inline void printAvailible() {
    #ifdef DEBUG
    for (int i = 0; i < W + 2; ++i) {
        for (int j = 0; j < H + 2; ++j) {
            eprint("%d", availible[i][j]);
        }
        eprint("\n");
    }
    #endif
}

void printMove(int x, int y) {
    eprint("my move %c %d\n", char(x - 1 + 'a'), y);
    print("move %c %d\n", char(x - 1 + 'a'), y);
}

int findMoveInDirection(int x0, int y0, int dx, int dy, int color,
    bool changeColor) {
    int x, y;
    for (x = x0 + dx, y = y0 + dy; map[x][y] == 1 - color;
        x += dx, y += dy
    ) {
    }
    if (!changeColor) {
        if (map[x][y] == EMPTY) {
            return max(abs(x - x0), abs(y - y0)) - 1;
        }
        return 0;
    }

    if (map[x][y] == color) {
        for (x = x0 + dx, y = y0 + dy; map[x][y] == 1 - color;
            x += dx, y += dy
        ) {
            if (color != myColor) {
                myStack.popElement(myStack.antiPos[x][y]);
            } else {
                myStack.push(x, y);
            }
            map[x][y] = color;
        }
    }
    return 0;
}

void recalcAvailible() {
    for (int x = 1; x <= W; ++x) {
        for (int y = 1; y <= H; ++y) {
            availible[x][y] = 0;
        }
    }

    for (int i = 0; i < myStack.ptr; ++i) {
        int x = myStack.xPos[i];
        int y = myStack.yPos[i];
        for (int j = 0; j < W; ++j) {
            int dx = directionsX[j];
            int dy = directionsY[j];

            int delta = findMoveInDirection(x, y, dx, dy, myColor, false);
            if (delta != 0) {
                availible[x + dx * (delta + 1)][y + dy * (delta + 1)] += delta;
            }
        }
    }

}

void makeMove(int x, int y, int color) {
    map[x][y] = color;
    for (int i = 0; i < W; ++i) {
        findMoveInDirection(x, y, directionsX[i], directionsY[i], color, true);
    }
    if (color == myColor) {
        myStack.push(x, y);
    }
    recalcAvailible();
}

void turn() {
    int maxi = 1;
    int maxj = 1;
    for (int i = 1; i <= W; ++i) {
        for (int j = 1; j <= W; ++j) {
            if (availible[i][j] >= availible[maxi][maxj]) {
                maxi = i;
                maxj = j;
            }
        }
    }
    makeMove(maxi, maxj, myColor);
    printMove(maxi, maxj);
}

void init() {
    for (int i = 0; i < W + 2; ++i) {
        map[0][i] = map[i][0] = map[W + 1][i] = map[i][W + 1] = OUTER;
    }
    for (int i = 1; i <= W; ++i) {
        for (int j = 1; j <= W; ++j) {
            map[i][j] = EMPTY;
        }
    }

    map[4][5] = WHITE;
    map[5][4] = WHITE;
    map[5][5] = BLACK;
    map[4][4] = BLACK;

    if (myColor == WHITE) {
        myStack.push(5, 4);
        myStack.push(4, 5);
    } else {
        myStack.push(5, 5);
        myStack.push(4, 4);
    }
    int lastDir = 0;
    for (int x = -1; x <= 1; ++x) {
        for (int y = -1; y <= 1; ++y) {
            if (x * x + y * y > 0) {
                directionsX[lastDir] = x;
                directionsY[lastDir] = y;
                ++lastDir;
            }
        }
    }
    recalcAvailible();
    printBoard();
    printAvailible();
}

int main() {
    char inputBuffer[MAX_BUFFER];
    scanf("init %s", inputBuffer);
    eprint("INITED\n");
    if (inputBuffer[0] == 'w') {
        myColor = 0;
    } else {
        myColor = 1;
    }
    init();
    while (true) {
        scanf("%s", inputBuffer);
        if (inputBuffer[0] == 't') {
            eprint("t\n");
            turn();
        } else if (inputBuffer[0] == 'm') {
            char c;
            int y;
            scanf(" %c %d", &c, &y);
            eprint("m\n");
            makeMove(c - 'a' + 1, y, 1 - myColor);
        } else {
            return 0;
        }
        printBoard();
        printAvailible();
    }
}
