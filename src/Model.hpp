#ifndef Model_hpp
#define Model_hpp

#include <map>
#include <string>
#include <vector>
class Alignment;
class Chunk;
class Parameter;
class Restaurant;
class StateSets;
class Tree;
class UserSettings;



class Model {

    public:
                                        Model(void) = delete;
                                        Model(Model& m) = delete;
                                        Model(UserSettings* s, int nss);
                                        Model(Alignment* aln, UserSettings* s);
                                       ~Model(void);
        Alignment*                      getAlignment(void) { return alignment; }
        std::vector<double>&            getBaseFrequencies(int subsetId);
        Chunk*                          getChunk(int idx) { return chunks[idx]; }
        int                             getDegreeTree(void);
        int                             getDegreeTreeLength(void);
        int                             getDegreeGammaShape(void);
        int                             getDegreeBaseFrequencies(void);
        int                             getDegreeRates(void);
        double                          getAlphaTree(void);
        double                          getAlphaTreeLength(void);
        double                          getAlphaShape(void);
        double                          getAlphaFrequencies(void);
        double                          getAlphaRates(void);
        std::string                     getPartitionTree(void);
        std::string                     getPartitionTreeLength(void);
        std::string                     getPartitionGammaShape(void);
        std::string                     getPartitionBaseFrequencies(void);
        std::string                     getPartitionRates(void);
        std::vector<double>&            getExchangeabilityRates(int subsetId);
        std::vector<double>&            getGammaCategoryRates(int subsetId);
        double                          getGammaShape(int subsetId);
        void                            getHeader(std::string& h);
        StateSets*                      getStateSetPtr(void) { return stateSets; }
        Tree*                           getTree(int id);
        double                          getTreeLength(void);
        double                          getTreeLength(int subsetId);
        std::vector<std::string>        getTaxonNames(void);
        void                            getParameterValues(std::string& v);
        double                          lnLikelihood(void);
        double                          lnLikelihood(int idx);
        Alignment*                      simulate(std::string fn, int ns);
        double                          update(void);
    
    private:
        char                            charCode(int x);
        void                            initializeDataChunks(UserSettings* s);
        void                            initializeParameters(UserSettings* s, std::vector<std::string> tn);
        void                            initializeProposalProbabilities(void);
        void                            initializeStateSets(UserSettings* s);
        Alignment*                      alignment;
        int                             numSubsets;
        StateSets*                      stateSets;
        std::vector<Chunk*>             chunks;
        std::vector<Parameter*>         parameters;
        Restaurant*                     treeRestaurant;
        Restaurant*                     treeLengthRestaurant;
        Restaurant*                     freqsRestaurant;
        Restaurant*                     ratesRestaurant;
        Restaurant*                     shapeRestaurant;
        std::map<Restaurant*,double>    proposalProbabilities;
};

#endif
