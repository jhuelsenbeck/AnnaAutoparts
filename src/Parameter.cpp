#include "Parameter.hpp"


bool Parameter::anyLessThanMin(std::vector<double>& f, double minVal) {

    for (int i=0; i<f.size(); i++)
        {
        if (f[i] < minVal)
            return true;
        }
    return false;
}
