#include <iomanip>
#include <iostream>
#include "UpdateInfo.hpp"



void UpdateInfo::accept(void) {
    if (attemptedUpdate == UpdateType::NO_MOVE)
        Msg::error("Tried to accept a move when there was none!");

    auto& moveStats = stats[attemptedUpdate];
    moveStats.numAccepted++;
    moveStats.numTries++;
        
    attemptedUpdate = UpdateType::NO_MOVE;
}

void UpdateInfo::reject(void) {
    if (attemptedUpdate == UpdateType::NO_MOVE)
        Msg::error("Tried to accept a move when there was none!");

    auto& moveStats = stats[attemptedUpdate];
    moveStats.numTries++;
        
    attemptedUpdate = UpdateType::NO_MOVE;
}

void UpdateInfo::resetStats(void) {
    for(auto& s : stats){
        s.numAccepted = 0;
        s.numTries = 0;
    }
}

void UpdateInfo::summary(void) {

    const auto& updateNames = UpdateNames();

    int paramWidth = 0;
    for (const auto& name : updateNames)
        paramWidth = std::max(paramWidth, (int)name.size());

    int triesWidth = 8;
    int accWidth = 8;
    int freqWidth  = 8;

    std::cout << "   * MCMC Acceptance Rates:\n";
    std::cout << "   * "
              << std::left << std::setw(paramWidth) << "Parameter"
              << " "
              << std::right << std::setw(triesWidth) << "Tries"
              << std::setw(accWidth)                 << "Accep"
              << std::setw(freqWidth)                << "Freq"
              << "\n";

    std::cout << "   * "
              << std::string(paramWidth, '-')
              << " "
              << std::string(triesWidth + accWidth + freqWidth, '-')
              << "\n";

    for (int i = 0; i < updateNames.size(); i++)
        {
        const auto& name = updateNames[i];
        auto& s = stats[i];

        std::cout << "   * "
                  << std::left  << std::setw(paramWidth) << name
                  << " "
                  << std::right << std::setw(triesWidth) << s.numTries
                  << std::setw(accWidth)                 << s.numAccepted
                  << std::fixed << std::setprecision(1)
                  << std::setw(freqWidth-1)              // Minus 1 for % sign
                  << (s.numTries > 0 ?
                      (double)s.numAccepted / s.numTries * 100.0 : 0.0)
                  << "%\n";
        }

    std::cout << std::endl;
}


void UpdateInfo::tune(){
    for(int i = 0; i < 9; i++)
        {
        auto& updateStatistic = stats[i];
        if(updateStatistic.numTries > 0)
            {
            double acceptRate = (double)updateStatistic.numAccepted/(double)updateStatistic.numTries;
            if(i == 5 || i == 6) {} // We skip local and TBR. These don't appear to tune well like this.
            else if(i > 0 && i <= 7)
                { // Handle moves with concentration type parameters
                if(acceptRate > 0.33) 
                    {
                    tunableParams[i] /= (1.0 + ((acceptRate-0.33)/0.67));
                    }
                else 
                    {
                    tunableParams[i] *= (2.0 - acceptRate/0.33);
                    }
                }
            else
                { // Handle moves with scale type parameters
                if (acceptRate > 0.33) 
                    {
                    tunableParams[i] *= (1.0 + ((acceptRate-0.33)/0.67));
                    }
                else 
                    {
                    tunableParams[i] /= (2.0 - acceptRate/0.33);
                    }
                }
            updateStatistic.numAccepted = 0;
            updateStatistic.numTries = 0;
            }
        }
}