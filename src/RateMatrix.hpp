#ifndef RateMatrix_hpp
#define RateMatrix_hpp

#include <vector>
#include "EigenSystem.hpp"



class RateMatrix {

    public:
                                                RateMatrix(void) = delete;
                                                RateMatrix(int na);
        void                                    updateRateMatrix(std::vector<double>& r, std::vector<double>& f);
    
    private:
        int                                     numMatrices;
        int                                     activeMatrix;
        std::vector<EigenSystem>                eigs;
};

#endif
