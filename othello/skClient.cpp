#define _CRT_SECURE_NO_WARNINGS

#include <cstdio>
#include <vector>
#include <algorithm>
#include <cstdarg>
#include <cassert>
#include <climits>
#include <iostream>
#include <ctime>

using std::vector;
using std::max;
using std::min;
using std::swap;
using std::abs;

clock_t start;

bool isTL() {
    return (((float)(clock() - start) / CLOCKS_PER_SEC) > 2.8);
}

const unsigned int MAXN = 65536;
const unsigned int MAXN_1 = 65535;
const unsigned int LOG_MAXN = 16;
const int MAX_DEPTH = 4;

const int WHITE = 0;
const int BLACK = 1;
const int EMPTY = 2;
const int OUTER = 4;

FILE *err;

#define stderr err

void makeMove(int, int, int, int);

const unsigned int MAX_BUFFER = 100;

const int W = 8, H = 8;

int map[W + 2][H + 2];
int availible[2][W + 2][H + 2];
int totalAvailible[2];
int totalScore[2];
int directionsX[W];
int directionsY[H];

void recalcAvailible(int color);

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

    void popElement(int x, int y) {
        popElement(antiPos[x][y]);
    }

    bool empty() {
        return ptr == 0;
    }
};

int myColor;

FastPositionsStack stack[2];

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
    for (int z = 0; z < 2; ++z) {
        eprint("----------------\n");
        for (int i = 0; i < W + 2; ++i) {
            for (int j = 0; j < H + 2; ++j) {
                eprint("%d ", availible[z][i][j]);
            }
            eprint("\n");
        }
    }
    #endif
}

void printMove(int x, int y) {
    eprint("my move %c %d\n", char(x - 1 + 'a'), y);
    print("move %c %d\n", char(x - 1 + 'a'), y);
}

void doBackUpAv(int backUpAvailible[2][W][H]) {
    for (int z = 0; z < 2; ++z) {
        for (int i = 0; i < W; ++i) {
            for (int j = 0; j < H; ++j) {
                backUpAvailible[z][i][j] = availible[z][i + 1][j + 1];
            }
        }
    }
}

void restoreBackUpAv(int backUpAvailible[2][W][H]) {
    for (int z = 0; z < 2; ++z) {
        for (int i = 0; i < W; ++i) {
            for (int j = 0; j < H; ++j) {
                 availible[z][i + 1][j + 1] = backUpAvailible[z][i][j];
                 totalAvailible[z] += bool(availible[z][i + 1][j + 1]);
            }
        }
    }
}

FastPositionsStack changedPositions;

int moves = 0;

int value[W][H] = {{13, 8, 11, 10, 10, 11, 8, 13},
                   {8, 8, 8, 8, 8, 8, 8, 8},
                   {11, 8, 8, 8, 8, 8, 8, 11},
                   {10, 8, 8, 8, 8, 8, 8, 10},
                   {10, 8, 8, 8, 8, 8, 8, 10},
                   {11, 8, 8, 8, 8, 8, 8, 11},
                   {8, 8, 8, 8, 8, 8, 8, 8},
                   {13, 8, 11, 10, 10, 11, 8, 13}};

void setCell(int x, int y, int color) {
    #ifdef DEBUG
        assert(map[x][y] == EMPTY);
    #endif
    map[x][y] = color;
    stack[color].push(x, y);
    totalScore[color] += value[x - 1][y - 1];
}

void clearCell(int x, int y) {
    #ifdef DEBUG
        assert(map[x][y] == WHITE || map[x][y] == BLACK);
    #endif
    stack[map[x][y]].popElement(x, y);
    totalScore[map[x][y]] -= value[x - 1][y - 1];
    map[x][y] = EMPTY;
}

void recolorCell(int x, int y) {
    #ifdef DEBUG
        assert(map[x][y] == WHITE || map[x][y] == BLACK);
    #endif
    totalScore[map[x][y]] -= value[x - 1][y - 1];
    totalScore[1 - map[x][y]] += value[x - 1][y - 1];
    stack[map[x][y]].popElement(x, y);
    stack[1 - map[x][y]].push(x, y);
    map[x][y] = 1 - map[x][y];
}


int minMaxFunction(int color) {
    if (!totalAvailible[0] && !totalAvailible[1]) {
        if (totalScore[color] > totalScore[1 - color]) {
            return INT_MAX / 4;
        } else if (totalScore[color] < totalScore[1 - color]) {
            return INT_MIN / 4;
        }
        return 0;
    }
    int a = 7, b = 13;
    if (moves > 40) {
        swap(a, b);
    }
    return 16 * b * (totalAvailible[color] - totalAvailible[1 - color]) + a * (totalScore[color] - totalScore[1 - color]);
}

int whereCurrentMaxX;
int whereCurrentMaxY;

int makeBestMove(int lastChangesPositions, int depth, int color) {
    if (depth == MAX_DEPTH || (!totalAvailible[0] && !totalAvailible[1])) {
        return minMaxFunction(color);
    }
    int currentMax = INT_MIN / 2;
    int backUp[2][W][H];
    doBackUpAv(backUp);

    int sign = 1;
    if (!totalAvailible[color]) {
        color = 1 - color;
        sign = -1;
    }

    for (int x = 1; x <= W; ++x) {
        for (int y = 1; y <= W; ++y) {
            if (availible[color][x][y]) {
                if (isTL() && depth > 1) {
                    return currentMax;
                }
                makeMove(x, y, color, 2);
                int value = -makeBestMove(changedPositions.ptr, depth + 1, 1 - color);
                if (value > currentMax) {
                    currentMax = value;
                    if (depth == 0) {
                        whereCurrentMaxX = x;
                        whereCurrentMaxY = y;
                    }
                }

                for (; changedPositions.ptr != lastChangesPositions;) {
                    recolorCell(changedPositions.topX(), changedPositions.topY());
                    changedPositions.pop();
                }
                clearCell(x, y);
                restoreBackUpAv(backUp);
            }
        }
    }
    eprint("currentMax: %d\n", currentMax);
    return sign * currentMax;
}


int findMoveInDirection(int x0, int y0, int dx, int dy, int color,
    int changeColor) {
    int x = x0 + dx, y = y0 + dy;
    for (; map[x][y] == 1 - color;
        x += dx, y += dy
    ) {
    }
    int result = max(abs(x - x0), abs(y - y0)) - 1;

    if (changeColor == 0) {
        if (map[x][y] == EMPTY) {
            return result;
        }
        return 0;
    }

    if (map[x][y] == color) {
        for (x = x0 + dx, y = y0 + dy; map[x][y] == 1 - color;
            x += dx, y += dy
        ) {
            if (changeColor == 2) {
                changedPositions.push(x, y);
            }
            recolorCell(x, y);
        }
        return result;
    }
    return 0;
}

void recalcAvailible(int color) {
    for (int x = 1; x <= W; ++x) {
        for (int y = 1; y <= H; ++y) {
            availible[color][x][y] = 0;
        }
    }
    totalAvailible[color] = 0;

    for (int i = 0; i < stack[color].ptr; ++i) {
        int x = stack[color].xPos[i];
        int y = stack[color].yPos[i];
        for (int j = 0; j < W; ++j) {
            int dx = directionsX[j];
            int dy = directionsY[j];

            int delta = findMoveInDirection(x, y, dx, dy, color, 0);
            if (delta != 0) {
                if (!availible[color][x + dx * (delta + 1)][y + dy * (delta + 1)]) {
                    ++totalAvailible[color];
                }
                availible[color][x + dx * (delta + 1)][y + dy * (delta + 1)] += delta;
            }
        }
    }
}

void makeMove(int x, int y, int color, int change=1) {
    setCell(x, y, color);
    for (int i = 0; i < W; ++i) {
        findMoveInDirection(x, y, directionsX[i], directionsY[i], color, change);
    }
    recalcAvailible(color);
    recalcAvailible(1 - color);
}

void turn() {
    makeBestMove(0, 0, myColor);
    makeMove(whereCurrentMaxX, whereCurrentMaxY, myColor);
    printMove(whereCurrentMaxX, whereCurrentMaxY);
    return;
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
    err = fopen("log.txt", "w");
    for (int i = 0; i < W + 2; ++i) {
        map[0][i] = map[i][0] = map[W + 1][i] = map[i][W + 1] = OUTER;
    }
    for (int i = 1; i <= W; ++i) {
        for (int j = 1; j <= W; ++j) {
            map[i][j] = EMPTY;
        }
    }

    setCell(4, 5, BLACK);
    setCell(5, 4, BLACK);
    setCell(4, 4, WHITE);
    setCell(5, 5, WHITE);

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
    recalcAvailible(WHITE);
    recalcAvailible(BLACK);
    printBoard();
    printAvailible();
}

int main() {
    char inputBuffer[MAX_BUFFER];
    scanf("init %s", inputBuffer);
    eprint("INITED\n");
    if (inputBuffer[0] == 'w') {
        myColor = WHITE;
    } else {
        myColor = BLACK;
    }
    init();
    FILE *log = fopen("log2.txt", "w");
    while (true) {
        scanf("%s", inputBuffer);
        if (inputBuffer[0] == 't') {
            eprint("t\n");
            fprint(log, "turn\n");
            start = clock();
            turn();
        } else if (inputBuffer[0] == 'm') {
            char c;
            int y;
            scanf(" %c %d", &c, &y);
            fprint(log, "move %c %d\n", c, y);
            eprint("m\n");
            makeMove(c - 'a' + 1, y, 1 - myColor);
        } else {
            return 0;
        }
        printBoard();
        printAvailible();
        ++moves;
    }
}
