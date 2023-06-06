#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include "Msg.hpp"
#include "UserSettings.hpp"

#undef DEBUG_MODE



UserSettings::UserSettings(int argc,  char* argv[]) {

    executablePath = argv[0];

    // set up a vector of arguments
    std::vector<std::string> settings;
#   if defined(DEBUG_MODE)
    settings.push_back("-i");
    settings.push_back("/Users/johnh/Desktop/AutoParts2_data/hummer.in");
    settings.push_back("-o");
    settings.push_back("/Users/johnh/Desktop/AutoParts2_data/hummer.out");
    settings.push_back("-n");
    settings.push_back("100000");
    settings.push_back("-p");
    settings.push_back("1000");
    settings.push_back("-b");
    settings.push_back("100");
    settings.push_back("-g");
    settings.push_back("4");
    settings.push_back("-eLength");
    settings.push_back("2.0");
    settings.push_back("-eShape");
    settings.push_back("2.0");
    settings.push_back("-ePi");
    settings.push_back("2.0");
    settings.push_back("-eTheta");
    settings.push_back("2.0");
    settings.push_back("-si");
    settings.push_back("/Users/johnh/Desktop/AutoParts2_data/sim.in");
    argc = (int)settings.size() + 1;
#   else
    for (int i=1; i<argc; i++)
        {
        std::string arg = argv[i];
        settings.push_back(arg);
        }
#   endif

    // check the number of arguments
    if (argc == 1)
        {
        usage();
        std::cout << "   * Executable path                                       = \"" << executablePath << "\"" << std::endl;
        Msg::error("Expecting command line arguments");
        }
    if (argc % 2 == 0)
        {
        usage();
        std::cout << "   * Executable path                                       = \"" << executablePath << "\"" << std::endl;
        Msg::error("Expecting an odd number of arguments");
        }

    // initialize settings variables
    inputFile                     = "";
    treeFile                      = "";
    outputFile                    = "";
    simFile                       = "";
    numMcmcCycles                 = 400000;
    burnIn                        = 0;
    numGammaCategories            = 1;
    treeLengthMean                = 1.0;
    treeLengthSD                  = 1.0;
    brlenLambda                   = 10.0;
    shapeLambda                   = 2.0;
    printFrequency                = 1000;
    sampleFrequency               = 100;
    isConcentrationParameterFixed = true;
	priorConcMean                 = 3.0;
	priorConcVariance             = 1.0;
    etAlpha                       = 2.0;
    etPi                          = 2.0;
    etTheta                       = 2.0;
    etLength                      = 2.0;
    tuningLocal                   = 200.0;
    tuningTreeLength              = log(4.0);
    tuningBaseFrequencies         = 300.0;
    tuningExchangabilityRates     = 300.0;
    tuningGammaShape              = log(4.0);
    tuningHeat                    = 0.1;
    tuningBrlen                   = 200.0;
    
    // interpret the arguments
    std::string currentArg = "";
    for (int i=0; i<settings.size(); i++)
        {
        if (currentArg == "")
            currentArg = settings[i];
        else
            {
            if (currentArg == "-i")
                inputFile = settings[i];
            else if (currentArg == "-t")
                treeFile = settings[i];
            else if (currentArg == "-o")
                outputFile = settings[i];
            else if (currentArg == "-si")
                simFile = settings[i];
            else if (currentArg == "-n")
                numMcmcCycles = stoi(settings[i]);
            else if (currentArg == "-p")
                printFrequency = stoi(settings[i]);
            else if (currentArg == "-s")
                sampleFrequency = stoi(settings[i]);
            else if (currentArg == "-b")
                burnIn = stoi(settings[i]);
            else if (currentArg == "-g")
                numGammaCategories = stoi(settings[i]);
            else if (currentArg == "-lenMean")
                treeLengthMean = stof(settings[i]);
            else if (currentArg == "-lenSD")
                treeLengthSD = stof(settings[i]);
            else if (currentArg == "-lenLam")
                brlenLambda = stof(settings[i]);
            else if (currentArg == "-e")
                shapeLambda = stof(settings[i]);
             else if (currentArg == "-c")
                isConcentrationParameterFixed = yesNo(settings[i]);
            else if (currentArg == "-m")
                priorConcMean = stof(settings[i]);
            else if (currentArg == "-k")
                priorConcVariance = stof(settings[i]);
            else if (currentArg == "-ePi")
                etPi = stof(settings[i]);
            else if (currentArg == "-eTheta")
                etTheta = stof(settings[i]);
            else if (currentArg == "-eShape")
                etAlpha = stof(settings[i]);
            else if (currentArg == "-eLength")
                etLength = stof(settings[i]);
            else if (currentArg == "-tLength")
                tuningTreeLength = stof(settings[i]);
            else if (currentArg == "-tLocal")
                tuningLocal = stof(settings[i]);
            else if (currentArg == "-tShape")
                tuningGammaShape = stof(settings[i]);
            else if (currentArg == "-tFreqs")
                tuningBaseFrequencies = stof(settings[i]);
            else if (currentArg == "-tRates")
                tuningExchangabilityRates = stof(settings[i]);
            else if (currentArg == "-tTBR")
                tuningHeat = stof(settings[i]);
            else if (currentArg == "-tBrlen")
                tuningBrlen = stof(settings[i]);
            else
                Msg::error("Could not interpret argument " + settings[i]);
            currentArg = "";
            }
        }
}

double UserSettings::getAlphaT(void) {

    return treeLengthMean * treeLengthMean / (treeLengthSD * treeLengthSD);
}

double UserSettings::getBetaT(void) {

    return treeLengthMean / (treeLengthSD * treeLengthSD);
}

void UserSettings::print(void) {

    std::cout << std::fixed << std::setprecision(3);
    
    std::cout << "   User Settings" << std::endl;
    std::cout << "   * Executable path                                       = \"" << executablePath << "\"" << std::endl;
    std::cout << "   * Input file path and name                              = \"" << inputFile << "\"" << std::endl;
    std::cout << "   * Tree file path and name                               = \"" << treeFile << "\"" << std::endl;
    std::cout << "   * Output file path and name                             = \"" << outputFile << "\"" << std::endl;
    
    std::cout << "   * Tree length mean                                      = " << treeLengthMean << std::endl;
    std::cout << "   * Tree length standard deviation                        = " << treeLengthSD << std::endl;
    std::cout << "   * Exponential parameter for gamma shape                 = " << shapeLambda << std::endl;
    std::cout << "   * Number of gamma categories                            = " << numGammaCategories << std::endl;
    
    std::cout << "   * Concentration parameter is a random variable          = " << ((isConcentrationParameterFixed == true) ? "yes" : "no") << std::endl;
    std::cout << "   * Concentration parameter mean (when r.v.)              = " << priorConcMean << std::endl;
    std::cout << "   * Concentration parameter variance (when r.v.)          = " << priorConcVariance << std::endl;
    std::cout << "   * Expected number of tables for tree length             = " << etLength << std::endl;
    std::cout << "   * Expected number of tables for gamma shape             = " << etAlpha << std::endl;
    std::cout << "   * Expected number of tables for base frequencies        = " << etPi << std::endl;
    std::cout << "   * Expected number of tables for exchangability rates    = " << etTheta << std::endl;
    
    std::cout << "   * Number of MCMC cycles                                 = " << numMcmcCycles << std::endl;
    std::cout << "   * Sample Frequency                                      = " << sampleFrequency << std::endl;
    std::cout << "   * Print Frequency                                       = " << printFrequency << std::endl;
    std::cout << "   * Burn in                                               = " << burnIn << std::endl;
    std::cout << "   * MCMC tuning parameter for the tree topology parameter = " << tuningLocal << std::endl;
    std::cout << "   * MCMC tuning parameter for the gamma shape parameter   = " << tuningGammaShape << std::endl;
    std::cout << "   * MCMC tuning parameter for the base frequencies        = " << tuningBaseFrequencies << std::endl;
    std::cout << "   * MCMC tuning parameter for the substitution rates      = " << tuningExchangabilityRates << std::endl;
    std::cout << "   * MCMC tuning parameter for the tree length parameter   = " << tuningTreeLength << std::endl;
    std::cout << "   * MCMC tuning parameter for the TBR heat parameter      = " << tuningHeat << std::endl;
    std::cout << "   * MCMC tuning parameter for the branch proportions      = " << tuningBrlen << std::endl;
    std::cout << std::endl;
}

void UserSettings::usage(void) {

    std::cout << "   Program options for data input/output:" << std::endl;
    std::cout << "   * -i       : Input file name" << std::endl;
    std::cout << "   * -t       : Tree file name (for constraining the analysis to a fixed tree)" << std::endl;
    std::cout << "   * -o       : Output file name" << std::endl;
    std::cout << std::endl;
    
    std::cout << "   Program options for model parameters:" << std::endl;
    std::cout << "   * -lenMean : Tree length mean" << std::endl;
    std::cout << "   * -lenSD   : Tree length standard deviation" << std::endl;
    std::cout << "   * -lenLam  : Tree length exponential parameter" << std::endl;
    std::cout << "   * -e       : Exponential parameter for shape parameter describing ASRV" << std::endl;
    std::cout << "   * -g       : Number of gamma categories" << std::endl;
    std::cout << std::endl;
    
    std::cout << "   Program options for DPP:" << std::endl;
    std::cout << "   * -c       : Concentration parameter is fixed (yes) or a random variable (no)" << std::endl;
    std::cout << "   * -k       : Prior mean of the number of categories when the concentration parameter is fixed" << std::endl;
    std::cout << "   * -m       : Prior mean of the concentration parameter when it is a random variable" << std::endl;
    std::cout << "   * -v       : Prior variance of the concentration parameter when it is a random variable" << std::endl;
    std::cout << "   * -eLength : Expected number of tables for tree length when the concentration parameter is fixed" << std::endl;
    std::cout << "   * -eShape  : Expected number of tables for gamma shape when the concentration parameter is fixed" << std::endl;
    std::cout << "   * -ePi     : Expected number of tables for base frequencies when the concentration parameter is fixed" << std::endl;
    std::cout << "   * -eTheta  : Expected number of tables for exchangability rates when the concentration parameter is fixed" << std::endl;
    std::cout << std::endl;
    std::cout << "   Program options for MCMC:" << std::endl;
    std::cout << "   * -n       : Number of MCMC cycles" << std::endl;
    std::cout << "   * -p       : Print-to-screen frequency" << std::endl;
    std::cout << "   * -s       : Sample-to-file frequency" << std::endl;
    std::cout << "   * -b       : Fraction of samples to discard (burn-in)" << std::endl;
    std::cout << "   * -tLocal  : MCMC tuning parameter for the tree topology parameter" << std::endl;
    std::cout << "   * -tTBR    : MCMC tuning parameter for heat parameter for biased TBR proposals" << std::endl;
    std::cout << "   * -tBrlen  : MCMC tuning parameter for the branch proportions update" << std::endl;
    std::cout << "   * -tShape  : MCMC tuning parameter for the gamma shape parameter" << std::endl;
    std::cout << "   * -tFreqs  : MCMC tuning parameter for the base frequencies" << std::endl;
    std::cout << "   * -tRates  : MCMC tuning parameter for the substitution rates" << std::endl;
    std::cout << "   * -tLength : MCMC tuning parameter for the tree length parameter" << std::endl;
}

bool UserSettings::yesNo(std::string str) {

    if ( (str[0] == 'Y' || str[0] == 'y') && (str[1] == 'E' || str[1] == 'e') && (str[2] == 'S' || str[2] == 's') )
        return true;
    else if ( (str[0] == 'N' || str[0] == 'n') && (str[1] == 'O' || str[1] == 'o') )
        return false;
    Msg::error("Unknown string " + str);
    return false;
}

