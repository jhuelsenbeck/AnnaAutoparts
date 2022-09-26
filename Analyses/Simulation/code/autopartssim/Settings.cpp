#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include "Settings.h"



Settings::Settings(std::string pf) {

    // open the file
	std::ifstream fileStrm( pf.c_str() );
	if (!fileStrm) 
		{
		std::cerr << "Cannot open file \"" + pf + "\"" << std::endl;
		exit(1);
		}

    // read the file line-by-line
	std::string linestring = "";
    std::string parmName = "";
    std::string val = "";
    std::vector<double> parmVector;
    std::vector< std::vector<double> > parmVals;
    double x = 0.0;
    bool readingParm = false, readingValue = false, readingVector = false, readingVectorTypeParm = false;
    int line = 0;
	while( getline(fileStrm, linestring).good() )
		{
        readingParm = true;
        val = "";
        for (int i=0; i<linestring.size(); i++)
            {
            char ch = linestring[i];

            if ( ch == ',' || ch == ')' || ch == ']' )
                {
                std::istringstream buf(val);
				buf >> x;
                parmVector.push_back( x );
                val = "";
                if ( (readingVector == true && ch == ')') || (readingVector == false && readingVectorTypeParm == false) )
                    {
                    /*std::cout << parmVector.size() << " -- ";
                    for (int j=0; j<parmVector.size(); j++)
                        std::cout << parmVector[j] << " ";
                    std::cout << std::endl;*/
                    parmVals.push_back( parmVector );
                    }
                }

            if (ch == ' ')
                {
                //std::cout << "space" << std::endl;
                }
            else if (ch == '[')
                {
                parmName = val;
                val = "";
                parmVector.clear();
                readingParm = false;
                readingValue = true;
                //std::cout << "left bracket" << std::endl;
                }
            else if (ch == ']')
                {
                readingVectorTypeParm = false;
                readingValue = false;
                parmVector.clear();
                //std::cout << "right bracket" << std::endl;
                }
            else if (ch == '(')
                {
                readingVectorTypeParm = true;
                parmVector.clear();
                readingVector = true;
                //std::cout << "left parantheses" << std::endl;
                }
            else if (ch == ')')
                {
                readingVector = false;
                parmVector.clear();
                //std::cout << "right parantheses" << std::endl;
                }
            else if (ch == ',')
                {
                if (readingVector == false)
                    parmVector.clear();
                //std::cout << "comma" << std::endl;
                }
            else if (isdigit(ch) == true || ch == '.')
                {
                val += ch;
                //std::cout << "digit (" << ch << ")" << std::endl;
                }
            else
                {
                val += ch;
                //std::cout << "something else (" << ch << ")" << std::endl;
                }
            
                
            }
            
        // print the parameters
#       if 0
        std::cout << "Parameter name = \"" << parmName << "\"" << std::endl;
        for (int i=0; i<parmVals.size(); i++)
            {
            for (int j=0; j<parmVals[i].size(); j++)
                std::cout << parmVals[i][j] << " ";
            std::cout << std::endl;
            }
#       endif

        // store the parameters
        if (parmName == "tree_length")
            {
            for (int i=0; i<parmVals.size(); i++)
                for (int j=0; j<parmVals[i].size(); j++)
                    branchLengthPrior.push_back( parmVals[i][j] );
            }
        else if (parmName == "num_taxa")
            {
            numTaxa = (int)(parmVals[0][0]);
            }
        else if (parmName == "num_reps")
            {
            numReplicates = (int)(parmVals[0][0]);
            }
        else if (parmName == "gamma_shape")
            {
            for (int i=0; i<parmVals.size(); i++)
                for (int j=0; j<parmVals[i].size(); j++)
                    alpha.push_back( parmVals[i][j] );
            }
        else if (parmName == "number_sites")
            {
            for (int i=0; i<parmVals.size(); i++)
                for (int j=0; j<parmVals[i].size(); j++)
                    numSites.push_back( (int)(parmVals[i][j]) );
            }
        else if (parmName == "base_frequencies")
            {
            for (int i=0; i<parmVals.size(); i++)
                {
                std::vector<double> temp;
                for (int j=0; j<parmVals[i].size(); j++)
                    temp.push_back( parmVals[i][j] );
                baseFreqs.push_back( temp );
                }
            }
        else if (parmName == "substitution_rates")
            {
            for (int i=0; i<parmVals.size(); i++)
                {
                std::vector<double> temp;
                for (int j=0; j<parmVals[i].size(); j++)
                    temp.push_back( parmVals[i][j] );
                substitutionRates.push_back( temp );
                }
            }
            
        for (int i=0; i<parmVals.size(); i++)
            parmVals[i].clear();
        parmVals.clear();
        line++;
        }

    // close the file
    fileStrm.close();
    
#   if 1
    std::cout << "Number of replicates: " << numReplicates << std::endl;

    std::cout << "Number of taxa:       " << numTaxa << std::endl;
    
    std::cout << "Branch length prior:  ";
    for (int i=0; i<branchLengthPrior.size(); i++)
        std::cout << branchLengthPrior[i] << " ";
    std::cout << std::endl;

    std::cout << "Gamma shape:          ";
    for (int i=0; i<alpha.size(); i++)
        std::cout << alpha[i] << " ";
    std::cout << std::endl;

    std::cout << "Number of sites:      ";
    for (int i=0; i<numSites.size(); i++)
        std::cout << numSites[i] << " ";
    std::cout << std::endl;
    
    std::cout << "Number of sites:      ";
    for (int i=0; i<baseFreqs.size(); i++)
        {
        std::cout << "( ";
        for (int j=0; j<baseFreqs[i].size(); j++)
            std::cout << baseFreqs[i][j] << " ";
        std::cout << ") ";
        }
    std::cout << std::endl;

    std::cout << "Substitution rates:   ";
    for (int i=0; i<substitutionRates.size(); i++)
        {
        std::cout << "( ";
        for (int j=0; j<substitutionRates[i].size(); j++)
            std::cout << substitutionRates[i][j] << " ";
        std::cout << ") ";
        }
    std::cout << std::endl;
#   endif

}