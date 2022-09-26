#include <cmath>
#include <iomanip>
#include <iostream>
#include "Alignment.hpp"
#include "ConditionalLikelihoods.hpp"

#undef CHECK_CONDLIKES



ConditionalLikelihoods::ConditionalLikelihoods(Alignment* aln, int nc) {

    // initialize some important variables
    numTaxa           = aln->getNumTaxa();
    numSites          = aln->getNumSites();
    numNodes          = 2 * numTaxa - 2;
    numStates         = 4;
    numRateCategories = nc;
    
    // allocate the conditional likelihoods and set up an array to position ourselves
    // in the raw vector of conditional likelihoods
    int clSize = numNodes * numSites * numStates * numRateCategories;
    for (int m=0; m<2; m++)
        {
        clsRaw[m] = new double[clSize];
        for (int i=0; i<clSize; i++)
            clsRaw[m][i] = 0.0;
        }
        
    // set up an array to help position ourselves in the conditional likelihoods
    for (int m=0; m<2; m++)
        {
        cls[m] = new double*[numNodes];
        double* p = &clsRaw[m][0];
        for (int i=0; i<numNodes; i++)
            {
            cls[m][i] = p;
            p += (numSites * numStates * numRateCategories);
            }
        }
            
    // initialize the tip conditional likelihoods
    for (int m=0; m<2; m++)
        {
        
        for (int i=0; i<numTaxa; i++)
            {
            double* p = cls[m][i];
            for (int j=0; j<numSites; j++)
                {
                int nucCode = aln->matrixEntry(i, j);
                int possibleNucs[4];
                aln->getPossibleNucs(nucCode, possibleNucs);
                for (int k=0; k<numRateCategories; k++)
                    {
                    for (int s=0; s<4; s++)
                        {
                        if (possibleNucs[s] == 1)
                            p[s] = 1.0;
                        else
                            p[s] = 0.0;
                        }
                    p += 4;
                    }

                }
            }
        
        }

#   if defined(CHECK_CONDLIKES)
    std::vector<std::string> errors = checkConditionalLikelihoods(aln, 0.00001);
    if (errors.size() > 0)
        {
        std::cout << errors.size() << " errors found in conditional likelihoods" << std::endl;
        exit(1);
        }
    else
        std::cout << "No errors found in conditional likelihoods" << std::endl;
#   endif
}

ConditionalLikelihoods::~ConditionalLikelihoods(void) {

    for (int m=0; m<2; m++)
        {
        delete [] clsRaw[m];
        delete [] cls[m];
        }
}

std::vector<std::string> ConditionalLikelihoods::checkConditionalLikelihoods(Alignment* aln, double tolerance) {

    std::vector<std::string> errorList;
    for (int space=0; space<2; space++)
        {
        for (int i=0; i<numTaxa; i++)
            {
            double* cl = getConditionalLikelihoods(space, i);
            for (int j=0; j<numSites; j++)
                {
                int nucCode = aln->matrixEntry(i, j);
                int possibleNucs[4];
                aln->getPossibleNucs(nucCode, possibleNucs);
                for (int k=0; k<numRateCategories; k++)
                    {
                    for (int s=0; s<numStates; s++)
                        {
                        double diff = fabs( (double)possibleNucs[s]-cl[s] );
                        if (diff > tolerance)
                            errorList.push_back("Error in conditional likelihoods at (" + std::to_string(i) + "," + std::to_string(j) + "," + std::to_string(k) + "," + std::to_string(s) + ")");
                        }
                    cl += numStates;
                    }

                }
            }
        }
    return errorList;
}

void ConditionalLikelihoods::print(void) {

    for (int c=0; c<numSites; c++)
        {
        std::cout << std::setw(4) << c << " -- ";
        std::cout << std::fixed << std::setprecision(0);
        for (int n=0; n<numTaxa; n++)
            {
            double* p = cls[0][n];
            p += c * numStates * numRateCategories;
            for (int k=0; k<numRateCategories; k++)
                {
                for (int i=0; i<numStates; i++)
                    std::cout << p[i];
                p += 4;
                }
            std::cout << " ";
            }
        std::cout << "| ";
        std::cout << std::endl;
        }
}
