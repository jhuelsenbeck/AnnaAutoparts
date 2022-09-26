#ifndef Parameter_hpp
#define Parameter_hpp

#include <string>
#include <vector>
class Model;
class Table;
class UserSettings;



class Parameter {

    public:
                                    Parameter(void) = delete;
                                    Parameter(Model* m, UserSettings* s, std::string n) { modelPtr = m; userSettings = s; name = n; }
        virtual                    ~Parameter(void) { }
        virtual void                accept(void) = 0;
        virtual std::string         getHeader(void) = 0;
        virtual std::string         getName(void) { return name; }
        virtual double              getProposalProbability(void) { return proposalProbability; }
        virtual Table*              getTable(void) { return myTable; }
        virtual bool                getUpdateModifiesEigens(void) { return updateModifiesEigens; }
        virtual std::string         getValuesAsString(int precision) = 0;
        virtual std::vector<double> getValue(void) = 0;
        virtual double              lnProbability(void) = 0;
        virtual Parameter*          newRandomInstance(void) = 0;
        virtual void                print(void) = 0;
        virtual void                reject(void) = 0;
        virtual void                setName(std::string s) { name = s; }
        virtual void                setProposalProbability(double x) { proposalProbability = x; }
        virtual void                setTable(Table* t) { myTable = t; }
        virtual void                setUpdateModifiesEigens(bool tf) { updateModifiesEigens = tf; }
        virtual std::string         type(void) = 0;
        virtual double              update(void) = 0;
    
    protected:
        bool                        anyLessThanMin(std::vector<double>& f, double minVal);
        std::string                 name;
        double                      proposalProbability;
        bool                        updateModifiesEigens;
        Model*                      modelPtr;
        Table*                      myTable;
        UserSettings*               userSettings;
};

#endif
