#ifndef Mcmc_hpp
#define Mcmc_hpp

#include <chrono>
#include <fstream>
#include <string>
class Model;
class UserSettings;

typedef std::chrono::high_resolution_clock::time_point Timer;



class Mcmc {

    public:
                        Mcmc(Model* m, UserSettings* s);
        void            run(void);
        void            burnin(void);
    
    private:
        void            closeOutputFiles(void);
        std::string     formattedTime(Timer& t1, Timer& t2);
        void            openOutputFiles(void);
        void            printToScreen(int n, double lnL, Timer& timePt, Timer& start);
        void            sample(int n, double lnL);
        std::ofstream   parmStrm;
        std::ofstream   treeStrm;
        int             numCycles;
        int             burnInCycles;
        int             printFrequency;
        int             sampleFrequency;
        int             tuningFrequency;
        Model*          model;
        UserSettings*   settings;
};

#endif
