#include <iomanip>
#include <iostream>
#include "Msg.hpp"
#include "UpdateInfo.hpp"



void UpdateInfo::accept(void) {

    if (attemptedUpdateName == "")
        Msg::error("Empty attempt string");

    std::map<std::string,UpdateStats>::iterator it = stats.find(attemptedUpdateName);
    if (it == stats.end())
        {
        UpdateStats v;
        v.numTries = 1;
        v.numAccepted = 1;
        stats.insert( std::make_pair(attemptedUpdateName,v) );
        }
    else
        {
        it->second.numTries++;
        it->second.numAccepted++;
        }
        
    attemptedUpdateName = "";
}

void UpdateInfo::reject(void) {

    if (attemptedUpdateName == "")
        Msg::error("Empty attempt string");

    std::map<std::string,UpdateStats>::iterator it = stats.find(attemptedUpdateName);
    if (it == stats.end())
        {
        UpdateStats v;
        v.numTries = 1;
        v.numAccepted = 0;
        stats.insert( std::make_pair(attemptedUpdateName,v) );
        }
    else
        {
        it->second.numTries++;
        }
        
    attemptedUpdateName = "";
}

void UpdateInfo::summary(void) {

    // get the length of the longest key
    int len = 0;
    for (std::map<std::string,UpdateStats>::iterator it = stats.begin(); it != stats.end(); it++)
        {
        if (it->first.length() > len)
            len = (int)it->first.length();
        }
        
    std::cout << "   MCMC Acceptance Rates:" << std::endl;
    std::cout << "   * Parameter";
    for (int i=0; i<len-9; i++)
        std::cout << " ";
    std::cout << "      Tries   Accep    Freq" << std::endl;
    std::cout << "   * ---------";
    for (int i=0; i<len-9; i++)
        std::cout << "-";
    std::cout << "---------------------------" << std::endl;
    for (std::map<std::string,UpdateStats>::iterator it = stats.begin(); it != stats.end(); it++)
        {
        std::cout << "   * " << it->first;
        for (int i=0; i<len-it->first.length(); i++)
            std::cout << " ";
        std::cout << "    " << std::setw(7) << it->second.numAccepted << " ";
        std::cout << std::setw(7) << it->second.numTries << " ";
        std::cout << std::fixed << std::setprecision(1) << std::setw(6) << ((double)it->second.numAccepted/it->second.numTries)*100.0 << "%";
        std::cout << std::endl;
        }
    std::cout << std::endl;
}
