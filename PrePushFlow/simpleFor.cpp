#include <cstdio>
#include "bits/stdc++.h"

using namespace std;

struct Edge {
    long long to;
    long long capacity;
    long long flow;
    long long reverse;
    int number;
    Edge() {}
    Edge(long long to, long long capacity, long long flow, long long reverse, int number) : 
    to(to), capacity(capacity), flow(flow), reverse(reverse), number(number) {}
};

vector <vector <Edge> > g;
vector <long long> goodEdge;
long long S, T;
long long n;

vector <long long> h;
vector <long long> ex;

long long extra = 0;

void relabel(long long v) {
    long long currentH = INT_MAX;
    for (int i = 0; i < g[v].size(); ++i) {
        const Edge &e = g[v][i];
        if (e.capacity - e.flow > 0)
            currentH = min(currentH, h[e.to]);
    }
    h[v] = currentH + 1;
}

void push(long long v, Edge &e) {
    long long pushValue = min(e.capacity - e.flow, ex[v]);
    ex[v] -= pushValue;
    ex[e.to] += pushValue;
    e.flow += pushValue;
    g[e.to][e.reverse].flow -= pushValue;
}

void discharge(long long u) {
    while (ex[u] > 0) {
        if (goodEdge[u] == g[u].size()) {
            relabel(u);
            goodEdge[u] = 0;
        }
        else {
            Edge &e = g[u][goodEdge[u]];
            if (e.capacity - e.flow > 0 && h[e.to] + 1 == h[u]) {
                push(u, e);
            }
            else {
                ++goodEdge[u];
            }
        }
    }
}
bool sf(const Edge &e1, const Edge &e2) {
    return e1.number < e2.number;
}

int main() {
    long long m;
    scanf("%lld %lld\n", &n, &m);
    S = 0;
    T = n - 1;
    g.resize(n);
    ex.assign(n, 0);
    h.assign(n, 0);
    goodEdge.assign(n, 0);
    for (long long i = 0; i < m; ++i)
    {
        long long a, b, c;
        scanf("%lld %lld %lld", &a, &b, &c);
        --a, --b;
        g[a].push_back(Edge(b, c, 0, (long long)g[b].size(), (int)i));
        g[b].push_back(Edge(a, 0, 0, (long long)g[a].size() - 1, (int)m));
    }
    for (int i = 0; i < g[S].size(); ++i) {
        Edge &e = g[S][i];
        e.flow = e.capacity;
        g[e.to][e.reverse].flow = -e.capacity;
        if (!ex[e.to] && e.capacity)
            ++extra;
        
        ex[e.to] += e.capacity;
    }
    for (long long i = 0; i < n; ++i)
        h[i] = 0;
    h[S] = n;
    while (extra) {
        for (long long i = 1; i + 1 < n; ++i)
            discharge(i);
        extra = 0;
        for (long long i = 1; i + 1 < n; ++i)
            if (ex[i] > 0) {
                ++extra;
//                 prlong longf("%d\n", i);
            }
//         prlong longf("---\n");
    }
    long long maxFlow = 0;
    for (int i = 0; i < g[S].size(); ++i) {
        Edge const &e = g[S][i];
        maxFlow += e.flow;
    }
    printf("%lld\n", maxFlow);
    vector <Edge> edges;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < g[i].size(); ++j) {
            edges.push_back(g[i][j]);
        }
    }
    sort(edges.begin(), edges.end(), sf);
    for (int i = 0; i < m; ++i) {
        printf("%lld\n", edges[i].flow);
    }
    return 0;
}
