#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include "Mcmc.hpp"
#include "Model.hpp"
#include "Msg.hpp"
#include "RandomVariable.hpp"
#include "Tree.hpp"
#include "UpdateInfo.hpp"
#include "UserSettings.hpp"



Mcmc::Mcmc(Model* m, UserSettings* s) {

    model = m;
    settings = s;
    
    numCycles = settings->getNumMcmcCycles();
    burnInCycles = settings->getBurnIn();
    tuningFrequency = settings->getTuningFrequency();
    printFrequency = settings->getPrintFrequency();
    sampleFrequency = settings->getSampleFrequency();
}

void Mcmc::closeOutputFiles(void) {

    parmStrm.close();
    treeStrm.close();
}

std::string Mcmc::formattedTime(Timer& t1, Timer& t2) {

    std::chrono::duration<double> durationSecs  = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1);
    int s = (int)durationSecs.count();
    int m = s / 60;
    int h = s / 3600;
        
    std::string tStr = "";
    if (h > 0)
        {
        tStr += std::to_string(h) + "h:";
        m -= h * 60;
        s -= h * 60 * 60;
        }
    if (m > 0 || (m == 0 && h > 0))
        {
        tStr += std::to_string(m) + "m:";
        s -= m * 60;
        }
    tStr += std::to_string(s) + "s";
    
    return tStr;
}

void Mcmc::openOutputFiles(void) {

    // open files for samples
    std::string outPath = settings->getOutputFile();
    std::string parmFileName = outPath + ".tsv";
    std::string treeFileName = outPath + ".tre";

    parmStrm.open( parmFileName.c_str(), std::ios::out );
    if (!parmStrm)
        Msg::error("Cannot open file \"" + parmFileName + "\"");
    treeStrm.open( treeFileName.c_str(), std::ios::out );
    if (!treeStrm)
        Msg::error("Cannot open file \"" + treeFileName + "\"");
}

void Mcmc::printToScreen(int n, double lnL, Timer& t2, Timer& t1) {

    std::cout << "   * " << std::setw(6) << n << " -- " << std::fixed << std::setprecision(2) << lnL << " -- ";
    std::cout << "tt(" << model->getDegreeTree() << ") ";
    std::cout << "tl(" << model->getDegreeTreeLength() << ") ";
    std::cout << "rv(" << model->getDegreeGammaShape() << ") ";
    std::cout << "bf(" << model->getDegreeBaseFrequencies() << ") ";
    std::cout << "er(" << model->getDegreeRates() << ") -- ";

    std::chrono::duration<double> durationSecs  = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1);
    double timePerCycle = (double)durationSecs.count() / n;
    if (timePerCycle == 0)
        timePerCycle = 1.0 / printFrequency;
    int s = (int)((numCycles - n) * timePerCycle);
    int m = s / 60;
    int h = s / 3600;

    if (h > 0)
        {
        std::cout << h << "h";
        std::cout << ":";
        m -= h * 60;
        s -= h * 60 * 60;
        }
    if (m > 0 || (m == 0 && h > 0))
        {
        if (m < 10)
            std::cout << "0";
        std::cout << m << "m";
        std::cout << ":";
        s -= m * 60;
        }
    if (s < 10)
        std::cout << "0";
    std::cout << s << "s";
    std::cout << " remaining  ";
    
    std::cout << std::endl;
}

void Mcmc::run(void) {

    std::cout << "   Markov Chain Monte Carlo Sampling:" << std::endl;
    // open files
    openOutputFiles();

    // run the chain
    auto start = std::chrono::high_resolution_clock::now();
    for (int n=1; n<=numCycles; n++)
        {
        // propose a new state for the chain
        double lnL = model->update();
        
        if (n % printFrequency == 0)
            {
            auto timePt = std::chrono::high_resolution_clock::now();
            printToScreen(n, lnL, timePt, start);
            }
                    
        if (n == 1 || n % sampleFrequency == 0)
            sample(n, lnL);
        }
    std::cout << std::endl;
    UpdateInfo::updateInfo().summary();
    
    // close file
    closeOutputFiles();
}

void Mcmc::burnin(void) {
    std::cout << "   Markov Chain Monte Carlo Burn-In:" << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    for (int n=1; n<=burnInCycles; n++)
        {
        // propose a new state for the chain
        double lnL = model->update();
        
        if(n % printFrequency == 0 || n == 1){
            std::cout << "   *    " << n << " -- " << lnL << std::endl;
        }

        if (n % tuningFrequency == 0)
            {
                std::cout << "\n   Tuning Parameters..." << std::endl;
                UpdateInfo::updateInfo().summary();
                UpdateInfo::updateInfo().tune();
            }
        }
    std::cout << std::endl;
    UpdateInfo::updateInfo().resetStats();
}

void Mcmc::sample(int n, double lnL) {

    // sample the tree
    Tree* t = model->getTree(0);
    if (n == 1)
        {
        treeStrm << "begin trees;" << std::endl;
        treeStrm << "   translate" << std::endl;
        std::vector<std::string> taxonNames = t->getTaxonNames();
        for (int i=0; i<taxonNames.size(); i++)
            {
            treeStrm << "   " << i+1 << " " << taxonNames[i];
            if (i == taxonNames.size()-1)
                treeStrm << ";" << std::endl;
            else
                treeStrm << "," << std::endl;
            }
        }
    treeStrm << "   tree t_" << n << " = " << t->getNewick() << std::endl;

    // sample the other parameters
    if (n == 1)
        {
        std::string parmHeader;
        model->getHeader(parmHeader);
        parmStrm << "Gen" << '\t' << "lnL" << '\t' << parmHeader << std::endl;
        }
    std::string parmStr;
    model->getParameterValues(parmStr);
    parmStrm << n << '\t' << lnL << '\t' << parmStr << std::endl;
}
