#include "ParameterSummary.hpp"
#include "ParameterSummaryPartition.hpp"
#include "ParameterSummaryReal.hpp"



bool ParameterSummary::isPartition(void) {

    ParameterSummaryPartition* p = dynamic_cast<ParameterSummaryPartition*>(this);
    if (p != NULL)
        return true;
    return false;
}

bool ParameterSummary::isReal(void) {

    ParameterSummaryReal* p = dynamic_cast<ParameterSummaryReal*>(this);
    if (p != NULL)
        return true;
    return false;
}
