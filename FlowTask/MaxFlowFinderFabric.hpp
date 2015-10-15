#ifndef _MAX_FLOW_FINDER_FABRIC_
#define _MAX_FLOW_FINDER_FABRIC_

#include "MaxFlowFinder.hpp"
#include "PrePushFlowSimpleFor.hpp"
#include "MalhotraKumarMaheshwari.hpp"

class MaxFlowFinderFabric {
public:
    enum EAlgorithms {
        PRE_PUSH_FLOW_SIMPLE_FOR,
        MALHOTRA_KUMAR_MAHESHWARI,
    };
    
    static MaxFlowFinder *getMaxFlowFinder(EAlgorithms algorithm) {
        switch(algorithm) {
            case PRE_PUSH_FLOW_SIMPLE_FOR:
                return new PrePushFlowSimpleFor();
            case MALHOTRA_KUMAR_MAHESHWARI:
                return new MalhotraKumarMaheshwari();
        }
    };
};

#endif