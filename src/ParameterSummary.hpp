#ifndef ParameterSummary_hpp
#define ParameterSummary_hpp

#include <string>
class Model;



class ParameterSummary {

    public:
        virtual void        addValue(double x) = 0;
        virtual void        addValue(std::string s) = 0;
        std::string         getName(void) { return name; }
        virtual int         getNumValues(void) = 0;
        virtual bool        inCi(Model* m) = 0;
        bool                isPartition(void);
        bool                isReal(void);
        virtual double      mse(Model* m) = 0;
        virtual void        print(void) = 0;
        virtual void        print(Model* m) = 0;
        void                setName(std::string s) { name = s; }
        virtual std::string summarize(void) = 0;
    
    protected:
        std::string         name;
};

#endif
