#include <fstream>
#include <iostream>
#include <vector>
#include "Alignment.h"
#include "MbRandom.h"
#include "Node.h"
#include "Settings.h"
#include "Tree.h"




Alignment::Alignment(Settings* sp, Tree* tp, MbRandom* rp) {

    ranPtr      = rp;
    settingsPtr = sp;
    treePtr     = tp;
    
    simulate();
}

Alignment::~Alignment(void) {

    delete [] matrix[0];
    delete [] matrix;
}

char Alignment::convertToNuc(int x) {

    if (x == 0)
        return 'A';
    else if (x == 1)
        return 'C';
    else if (x == 2)
        return 'G';
    else 
        return 'T';
}

void Alignment::print(void) {

    for (int i=0; i<settingsPtr->getNumTaxa(); i++)
        {
        std::cout << "Taxon_" << i+1 << "   ";
        for (int j=0; j<totalLength; j++)
            std::cout << convertToNuc(matrix[i][j]);
        std::cout << std::endl;
        }
}

void Alignment::printToFile(std::string fn) {

    // open the output files
    std::string apFileName = fn + "_ap.in";
    std::string nexFileName = fn + ".nex";
    std::ofstream apStrm, nexStrm;
    apStrm.open( apFileName.c_str(), std::ios::out );
    if ( !apStrm )
        {
        std::cerr << "ERROR: Problem opening output file" << std::endl;
        exit(1);
        }
    nexStrm.open( nexFileName.c_str(), std::ios::out );
    if ( !nexStrm )
        {
        std::cerr << "ERROR: Problem opening output file" << std::endl;
        exit(1);
        }

    // output the autoparts file
    apStrm << settingsPtr->getNumTaxa() << " " << totalLength << std::endl;
    for (int i=0; i<settingsPtr->getNumTaxa(); i++)
        {
        apStrm << "Taxon_" << i+1 << "   ";
        if (i+1 < 10)
            apStrm << " ";
        for (int j=0; j<totalLength; j++)
            apStrm << convertToNuc(matrix[i][j]);
        apStrm << std::endl;
        }
    std::vector<int> numSites = settingsPtr->getNumSites();
    int begSite = 1, endSite = 0;
    for (int i=0; i<numSites.size(); i++)
        {
        endSite += numSites[i];
        apStrm << "charset " << begSite << "-" << endSite << ";" << std::endl;
        begSite = endSite + 1;
        }
   
    // output the nexus file for the maximally partitionned analysis
    nexStrm << "#NEXUS" << std::endl << std::endl;
    nexStrm << "begin data;" << std::endl;
    nexStrm << "   dimensions ntax=" << settingsPtr->getNumTaxa() << " nchar=" << totalLength << ";" << std::endl;
    nexStrm << "   format datatype=dna;" << std::endl;
    nexStrm << "   matrix" << std::endl;
    for (int i=0; i<settingsPtr->getNumTaxa(); i++)
        {
        nexStrm << "   Taxon_" << i+1 << "   ";
        if (i+1 < 10)
            nexStrm << " ";
        for (int j=0; j<totalLength; j++)
            nexStrm << convertToNuc(matrix[i][j]);
        nexStrm << std::endl;
        }
    nexStrm << "   ;" << std::endl;
    nexStrm << "end;" << std::endl << std::endl;
    
    nexStrm << "begin mrbayes;" << std::endl;
	nexStrm << "   set autoclose=yes;" << std::endl;
	nexStrm << "   log start filename=mb_sim2_Kmax_log.txt;" << std::endl;			// we'll want to increment the log-file name here
	begSite = 1;
    endSite = 0;
    for (int i=0; i<numSites.size(); i++)
        {
        endSite += numSites[i];
        nexStrm << "   charset part_" << i+1 << " = " << begSite << "-" << endSite << ";" << std::endl;
        begSite = endSite + 1;
        }
    nexStrm << "   partition my_part = " << numSites.size() << ":";
    for (int i=0; i<numSites.size(); i++)
        {
        nexStrm << "part_" << i+1;
        if ( i + 1 != numSites.size() )
            nexStrm << ",";
        }
    nexStrm << ";" << std::endl;
    nexStrm << "   set partition=my_part;" << std::endl;
    nexStrm << "   lset applyto=(all) nst=6 rates=gamma;" << std::endl;
    nexStrm << "   prset revmatpr=dirichlet(1,1,1,1,1,1) statefreqpr=dirichlet(1,1,1,1) shapepr=uniform(0.1,50);" << std::endl;
    nexStrm << "   unlink statefreq=(all) revmat=(all) shape=(all);" << std::endl;
	nexStrm << "   prset applyto=(all) ratepr=variable;" << std::endl;
	nexStrm << "   mcmc ngen=2000000 printfreq=1000 samplefreq=1000" << std::endl;
	nexStrm << "   nchains=4 savebrlens=yes filename=mb_sim2_Kmax;" << std::endl;	// we'll want to increment the log-file name here
	nexStrm << "   sumt filename=mb_sim2_Kmax burnin=500;" << std::endl;			// we'll want to increment the log-file name here
	nexStrm << "   sump filename=mb_sim2_Kmax burnin=500;" << std::endl;			// we'll want to increment the log-file name here
	nexStrm << "end;" << std::endl;

    // close the file streams
    apStrm.close();
    nexStrm.close();
}

void Alignment::setRateMatrix(double q[4][4], int whichPart) {

    // get the parameters
    std::vector< std::vector<double> > f = settingsPtr->getBaseFreqs();
    std::vector< std::vector<double> > r = settingsPtr->getSubstitutionRates();
    
    // set the off-diagonal parts of the rate matrix
    for (int i=0, k=0; i<4; i++)
        {
        for (int j=i+1; j<4; j++)
            {
            q[i][j] = r[whichPart][k] * f[whichPart][j];
            q[j][i] = r[whichPart][k] * f[whichPart][i];
            k++;
            }
        }
    
    // set the diagonals
    double averageRate = 0.0;
    for (int i=0; i<4; i++)
        {
        double sum = 0.0;
        for (int j=0; j<4; j++)
            {
            if (i != j)
                sum += q[i][j];
            }
        q[i][i] = -sum;
        averageRate += f[whichPart][i] * sum;
        }
        
    // rescale the rate matrix so the average rate is one
    double factor = 1.0 / averageRate;
    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            q[i][j] *= factor;
}

void Alignment::simulate(void) {

    // get the size of the alignment
    std::vector<int> numSites = settingsPtr->getNumSites();
    totalLength = 0;
    for (int i=0; i<numSites.size(); i++)
        totalLength += numSites[i];
    
    // allocate the matrix
    matrix = new int*[settingsPtr->getNumTaxa()];
    matrix[0] = new int[settingsPtr->getNumTaxa() * totalLength];
    for (int i=1; i<settingsPtr->getNumTaxa(); i++)
        matrix[i] = matrix[i-1] + totalLength;
    for (int i=0; i<settingsPtr->getNumTaxa(); i++)
        for (int j=0; j<totalLength; j++)
            matrix[i][j] = 0;
        
    // get some important parameters for the simulation
    std::vector<Node*> downPass = treePtr->getTraversalSequence();
    std::vector<double> brlen = settingsPtr->getBranchLengthPrior();
    std::vector<double> alpha = settingsPtr->getAlpha();
    std::vector< std::vector<double> > f = settingsPtr->getBaseFreqs();
    
    // simulate the data on a partition-by-partition basis
    int siteNum = 0;
    for (int part=0; part<numSites.size(); part++)
        {
        double q[4][4];
        setRateMatrix(q, part);
        double treeLength = brlen[part];
        for (int site=0; site<numSites[part]; site++)
            {
            double r = ranPtr->gammaRv(alpha[part], alpha[part]);
            for (std::vector<Node*>::reverse_iterator p = downPass.rbegin(); p != downPass.rend(); p++)
                {
                if ( (*p) == treePtr->getRoot() )
                    {
                    // draw from the stationary distribution
                    double u = ranPtr->uniformRv();
                    double sum = 0.0;
                    for (int i=0; i<4; i++)
                        {
                        sum += f[part][i];
                        if (u < sum)
                            {
                            (*p)->setNuc( i );
                            break;
                            }
                        }
                    }
                else 
                    {
                    // simulate along a branch of the tree
                    double v = (*p)->getProportion() * treeLength * r;
                    int curState = (*p)->getAnc()->getNuc();
                    double curBranchPos = 0.0;
                    while (curBranchPos < v)
                        {
                        double rate = -q[curState][curState];
                        curBranchPos += ranPtr->exponentialRv(rate);
                        if (curBranchPos < v)
                            {
                            double u = ranPtr->uniformRv();
                            double sum = 0.0;
                            for (int j=0; j<4; j++)
                                {
                                if (curState != j)
                                    {
                                    sum += q[curState][j] / rate;
                                    if (u < sum)
                                        {
                                        curState = j;
                                        break;
                                        }
                                    }
                                }
                            }
                        }
                    (*p)->setNuc(curState);
                    }
                if ( (*p)->getLft() == NULL || (*p)->getRht() == NULL || (*p)->getAnc() == NULL)
                    matrix[(*p)->getIndex()][siteNum] = (*p)->getNuc();
                }
            siteNum++;
            }
        }
}