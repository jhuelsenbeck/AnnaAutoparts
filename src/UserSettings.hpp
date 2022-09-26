#ifndef UserSettings_hpp
#define UserSettings_hpp

#include <string>



class UserSettings {

    public:
                        UserSettings(void) = delete;
                        UserSettings(int argc,  char* argv[]);
        double          getAlphaT(void);
        double          getBetaT(void);
        int             getBurnIn(void) { return burnIn; }
        double          getGammaShapeLambda(void) { return shapeLambda; }
        double          getExpectedNumberTreeLengthTables(void) { return etLength; }
        double          getExpectedNumberAlphaTables(void) { return etAlpha; }
        double          getExpectedNumberPiTables(void) { return etPi; }
        double          getExpectedNumberThetaTables(void) { return etTheta; }
        bool            getIsConcentrationParameterFixed(void) { return isConcentrationParameterFixed; }
        double          getPriorConcMean(void) { return priorConcMean; }
        double          getPriorConcVariance(void) { return priorConcVariance; }
        std::string     getInputFile(void) { return inputFile; }
        std::string     getOutputFile(void) { return outputFile; }
        int             getNumGammaCategories(void) { return numGammaCategories; }
        int             getNumMcmcCycles(void) { return numMcmcCycles; }
        int             getPrintFrequency(void) { return printFrequency; }
        int             getSampleFrequency(void) { return sampleFrequency; }
        std::string     getSimFile(void) { return simFile; }
        std::string     getTreeFile(void) { return treeFile; }
        double          getTuningBrlen(void) { return tuningBrlen; }
        double          getTuningLocal(void) { return tuningLocal; }
        double          getTuningTreeLength(void) { return tuningTreeLength; }
        double          getTuningBaseFrequencies(void) { return tuningBaseFrequencies; }
        double          getTuningExchangabilityRates(void) { return tuningExchangabilityRates; }
        double          getTuningGammaShape(void) { return tuningGammaShape; }
        double          getTuningHeat(void) { return tuningHeat; }
        void            print(void);
        void            setInputFile(std::string s) { inputFile = s; }
        void            setNumMcmcCycles(int x) { numMcmcCycles = x; }
        void            setOutputFile(std::string s) { outputFile = s; }
        void            usage(void);
        
    private:
        bool            yesNo(std::string str);
        std::string     inputFile;
        std::string     treeFile;
        std::string     outputFile;
        std::string     simFile;
        int             numMcmcCycles;
        int             printFrequency;
        int             sampleFrequency;
        int             burnIn;
        int             numGammaCategories;
        double          treeLengthMean;
        double          treeLengthSD;
        double          shapeLambda;
        bool            isConcentrationParameterFixed;
        double          priorConcMean;
        double          priorConcVariance;
        double          etLength;
        double          etAlpha;
        double          etPi;
        double          etTheta;
        double          tuningLocal;
        double          tuningHeat;
        double          tuningBrlen;
        double          tuningTreeLength;
        double          tuningBaseFrequencies;
        double          tuningExchangabilityRates;
        double          tuningGammaShape;
};

#endif
