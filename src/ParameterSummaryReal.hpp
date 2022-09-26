#ifndef ParameterSummaryReal_hpp
#define ParameterSummaryReal_hpp

#include <string>
#include <vector>
#include "ParameterSummary.hpp"


class ParameterSummaryReal : public ParameterSummary {

    public:
        void                        addValue(double x) { values.push_back(x); }
        void                        addValue(std::string s) { }
        int                         getNumValues(void) { return (int)values.size(); }
        bool                        inCi(Model* m);
        double                      mse(Model* m);
        void                        print(void);
        void                        print(Model* m);
        std::string                 summarize(void);
    
    private:
        std::pair<double,double>    calculateCredibleInterval(double cv);
        double                      calculateMean(void);
        double                      calculateVariance(void);
        double                      getTrueValue(Model* m);
        std::vector<double>         values;
};

#endif
