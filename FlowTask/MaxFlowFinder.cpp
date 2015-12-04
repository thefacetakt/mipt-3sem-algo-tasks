#include "MaxFlowDescription.hpp"
#include "Net.hpp"
#include "MaxFlowFinder.hpp"
#include <vector>
#include <climits>


MaxFlowDescription MaxFlowFinder::run(unsigned int verticesCount, const std::vector<DirectedEdgeWithStart> &net, unsigned int source, unsigned int sink) {
    if (source == sink) {
        MaxFlowDescription result;
        result.flowValue = ULLONG_MAX;
        return result;
    }
    net_ = Net(verticesCount, net, source, sink);
    init(verticesCount);
    findMaxFlow();
    MaxFlowDescription result = net_.returnFlowDescription();
    net_.cleanUp();
    cleanUp();
    return result;
}