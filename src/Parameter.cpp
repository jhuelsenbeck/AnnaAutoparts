#include <iostream>
#include "Parameter.hpp"



Parameter::Parameter(const Parameter& parm) {

    std::cout << "Parameter copy constructor" << std::endl;
    
    name = parm.name;
    proposalProbability = parm.proposalProbability;
    updateModifiesEigens = parm.updateModifiesEigens;
    modelPtr = parm.modelPtr;
    myTable = parm.myTable;
    userSettings = parm.userSettings;
}

bool Parameter::anyLessThanMin(std::vector<double>& f, double minVal) {

    for (int i=0; i<f.size(); i++)
        {
        if (f[i] < minVal)
            return true;
        }
    return false;
}
