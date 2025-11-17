#ifndef ParameterTreeLength_hpp
#define ParameterTreeLength_hpp

#include <vector>
#include "Parameter.hpp"
class Model;
class UserSettings;



class ParameterTreeLength : public Parameter {

    public:
                                ParameterTreeLength(void) = delete;
                                ParameterTreeLength(const ParameterTreeLength& parm);
                                ParameterTreeLength(Model* m, UserSettings* s, double a, double b);
        void                    accept(void);
        std::string             getHeader(void);
        std::vector<double>     getValue(void);
        std::string             getValuesAsString(int precision);
        double                  lnProbability(void);
        Parameter*              newRandomInstance(void);
        void                    print(void);
        void                    reject(void);
        std::string             type(void);
        double                  update(void);

    private:
        double                  alphaT;
        double                  betaT;
        double                  length[2];
};

#endif
