#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include "Alignment.hpp"
#include "Mcmc.hpp"
#include "Model.hpp"
#include "Msg.hpp"
#include "ParameterSummary.hpp"
#include "ParameterSummaryPartition.hpp"
#include "ParameterSummaryReal.hpp"
#include "UserSettings.hpp"

#define POSTERIOR_ANALYSIS

void printHeader(void);
void readTsvFile(std::string fn, int burn, std::vector<ParameterSummary*>& parms);



#if defined(POSTERIOR_ANALYSIS)



void summarize(std::string fn, int burn, std::vector<ParameterSummary*>& parms);



// posterior analysis of a data file
int main(int argc, char* argv[]) {

    // print header
    printHeader();

    // get the user settings
    UserSettings settings(argc, argv);
    settings.print();

    // read in the data
    Alignment data(settings.getInputFile());
    //data.print();
        
    // set up the phylogenetic model
    Model model(&data, &settings);
    
    // perform the MCMC analysis
    Mcmc mcmc(&model, &settings);
    mcmc.run();

    // summarize the results
    std::vector<ParameterSummary*> parms;
    summarize(settings.getOutputFile()+".tsv", settings.getBurnIn(), parms);

    return 0;
}

void summarize(std::string fn, int burn, std::vector<ParameterSummary*>& parms) {

    readTsvFile(fn, burn, parms);
    
    std::cout << "   Results Summary:" << std::endl;
    for (int i=0; i<parms.size(); i++)
        parms[i]->print();
}



#else



struct Results {

    int n;
    std::vector<std::string> name;
    std::vector<double> mse;
    std::vector<double> ci;
};

void summarize(std::string fn, int burn, std::vector<ParameterSummary*>& parms, Model* trueModel, Results& res);
void tabulateResults(std::vector<ParameterSummary*>& parms, Model* trueModel, Results& res);



// analysis of simulated data
int main(int argc, char* argv[]) {

    // print header
    printHeader();

    // get the user settings
    UserSettings settings(argc, argv);
    settings.print();
    
    int numTaxa = 10;
    int numReplicates = 100;
    int numPartitions = 5;
    int numSitesPerPartition = 100;
    Results results;
    results.n = 0;
    for (int i=0; i<numReplicates; i++)
        {
        // set up the phylogenetic model
        Model simModel(&settings, numTaxa, numPartitions);
        
        // get the simulated data
        Alignment* data = simModel.simulate(settings.getSimFile(), numSitesPerPartition);
                
        // set up the phylogenetic model
        //Model model(data, &settings);
        Model model(data, &settings, &simModel);
        
        // perform the MCMC analysis
        Mcmc mcmc(&model, &settings);
        mcmc.run();
        
        // compare to the true results
        std::vector<ParameterSummary*> parms;
        summarize(settings.getOutputFile()+".tsv", settings.getBurnIn(), parms, &simModel, results);
        
        delete data;
        }
        
    // summarize results over simulations
    for (int i=0; i<results.name.size(); i++)
        std::cout << results.name[i] << " -- " << results.mse[i]/numReplicates << " " << results.ci[i]/numReplicates << std::endl;
    
    return 0;
}

void summarize(std::string fn, int burn, std::vector<ParameterSummary*>& parms, Model* trueModel, Results& res) {

    readTsvFile(fn, burn, parms);
    tabulateResults(parms, trueModel, res);

    std::cout << "   Results Summary:" << std::endl;
    for (int i=0; i<parms.size(); i++)
        parms[i]->print(trueModel);
}

void tabulateResults(std::vector<ParameterSummary*>& parms, Model* trueModel, Results& res) {

    for (int i=0; i<parms.size(); i++)
        {
        double mse = parms[i]->mse(trueModel);
        bool inCi = parms[i]->inCi(trueModel);
        if (res.n == 0)
            {
            res.name.push_back(parms[i]->getName());
            res.mse.push_back(mse);
            res.ci.push_back((double)inCi);
            }
        else
            {
            res.mse[i] += mse;
            res.ci[i] += (double)inCi;
            }
        }
    if (res.n == 0)
        res.n = 1;
    else
        res.n++;
}

#endif

void printHeader(void) {

    std::cout << "   AutoParts v1.0" << std::endl;
    std::cout << "   * John P. Huelsenbeck(1), Brian Moore(2), and Anna Chriss(1)" << std::endl;
    std::cout << "   * (1) University of California, Berkeley" << std::endl;
    std::cout << "   * (2) University of California, Davis" << std::endl << std::endl;
}

void readTsvFile(std::string fn, int burn, std::vector<ParameterSummary*>& parms) {

    // open the file
    std::ifstream strm(fn.c_str());
    if (!strm)
        Msg::error("Cannot open file \"" + fn + "\"");

    std::string lineString = "";
    int line = 0;
    std::vector<std::string> header;
    std::vector<bool> isPartition;
    while( getline(strm, lineString).good() )
        {
        std::istringstream linestream(lineString);
        //std::cout << lineString << std::endl;
        int ch;
        std::string word = "";
        int wordNum = 0;
        std::string cmdString = "";
        do
            {
            word = "";
            linestream >> word;
            wordNum++;
            if (line == 0)
                {
                if (word != "")
                    {
                    bool isRgf = false;
                    if (word.at(0) == 'R' && word.at(1) == 'G' && word.at(2) == 'F')
                        isRgf = true;
                    isPartition.push_back(isRgf);

                    header.push_back(word);
                    if (isRgf == false)
                        {
                        ParameterSummaryReal* newParm = new ParameterSummaryReal;
                        parms.push_back(newParm);
                        newParm->setName(word);
                        }
                    else
                        {
                        ParameterSummaryPartition* newParm = new ParameterSummaryPartition;
                        parms.push_back(newParm);
                        newParm->setName(word);
                        }
                    }
                }
            else
                {
                
                if (line > burn + 1 && word != "")
                    {
                    if (isPartition[wordNum-1] == false)
                        {
                        double x;
                        std::istringstream buf(word);
                        buf >> x;
                        parms[wordNum-1]->addValue(x);
                        }
                    else
                        {
                        parms[wordNum-1]->addValue(word);
                        }
                    }
                    
                }
            } while ( (ch=linestream.get()) != EOF );
            
        line++;
        }
    
    // close the file
    strm.close();
}


