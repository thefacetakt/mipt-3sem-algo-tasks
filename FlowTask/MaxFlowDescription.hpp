#ifndef _MAX_FLOW_DESCRIPTION_
#define _MAX_FLOW_DESCRIPTION_

#include "Edges.hpp"
#include <vector>

class MaxFlowDescription {
public:
    unsigned long long flowValue;
    std::vector<FlowEdge> description;
};

#endif