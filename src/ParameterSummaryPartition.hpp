#ifndef ParameterSummaryPartition_hpp
#define ParameterSummaryPartition_hpp

#include <map>
#include <string>
#include <vector>
#include "ParameterSummary.hpp"
class Partition;



class ParameterSummaryPartition : public ParameterSummary {

    public:
                                    ParameterSummaryPartition(void);
                                   ~ParameterSummaryPartition(void);
        void                        addValue(double x) { }
        void                        addValue(std::string s);
        int                         getNumValues(void) { return 0; }
        bool                        inCi(Model* m);
        double                      mse(Model* m);
        void                        print(void);
        void                        print(Model* m);
        std::string                 summarize(void);
    
    private:
        void                        calculateMeanPartition(void);
        std::string                 getTruePartition(Model* m);
        void                        pairProbabilities(void);
        int                         numSamples;
        std::map<std::string,int>   partitionsStringMap;
        std::map<Partition*,int>    partitionsMap;
};

#endif
