#include "MaxFlowFinder.hpp"
#include "MaxFlowFinderFabric.hpp"
#include "Edges.hpp"
#include <cstdio>
#include <vector>

int main() {
    unsigned int n, m;
    scanf("%u %u", &n, &m);
    std::vector <DirectedEdgeWithStart> graph(m);
    for (unsigned int i = 0; i < m; ++i) {
        unsigned int from, to, capacity;
        scanf("%u %u %u", &from, &to, &capacity);
        --from, --to;
        graph[i] = DirectedEdgeWithStart(from, to, capacity);
    }
    MaxFlowFinder *maxFlowFinder = MaxFlowFinderFabric::getMaxFlowFinder(MaxFlowFinderFabric::PRE_PUSH_FLOW_SIMPLE_FOR);
    MaxFlowDescription description = maxFlowFinder->findMaxFlow(n, graph, 0, n - 1);
    printf("%llu\n",  description.flowValue);
    delete maxFlowFinder;
    return 0;
}