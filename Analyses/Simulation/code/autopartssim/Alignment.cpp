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

void Alignment::printToFile(std::string fn, int rep) {

	// get a string with the replicate number
	char cStr[20];
	sprintf(cStr, "%d", rep+1);
	std::string repStr = cStr;
	
    // open the output files
    std::string apFileName = fn + "_ap_" + repStr + ".in";
    std::string nexFile1Name = fn + "_mb_max_" + repStr + ".nex";
    std::string nexFile2Name = fn + "_mb_tru_" + repStr + ".nex";
    std::ofstream apStrm, nexStrm1, nexStrm2;
    apStrm.open( apFileName.c_str(), std::ios::out );
    if ( !apStrm )
        {
        std::cerr << "ERROR: Problem opening ap output file" << std::endl;
        exit(1);
        }
    nexStrm1.open( nexFile1Name.c_str(), std::ios::out );
    if ( !nexStrm1 )
        {
        std::cerr << "ERROR: Problem opening nex max output file" << std::endl;
        exit(1);
        }
    nexStrm2.open( nexFile2Name.c_str(), std::ios::out );
    if ( !nexStrm2 )
	{
        std::cerr << "ERROR: Problem opening nex tru output file" << std::endl;
        exit(1);
	}
	
    // output the phylip file for analysis under the DPP model
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
   
    // output the nexus file analysis under the maximally partitionned model
    nexStrm1 << "#NEXUS" << std::endl << std::endl;
    nexStrm1 << "begin data;" << std::endl;
    nexStrm1 << "   dimensions ntax=" << settingsPtr->getNumTaxa() << " nchar=" << totalLength << ";" << std::endl;
    nexStrm1 << "   format datatype=dna;" << std::endl;
    nexStrm1 << "   matrix" << std::endl;
    for (int i=0; i<settingsPtr->getNumTaxa(); i++)
        {
        nexStrm1 << "   Taxon_" << i+1 << "   ";
        if (i+1 < 10)
            nexStrm1 << " ";
        for (int j=0; j<totalLength; j++)
            nexStrm1 << convertToNuc(matrix[i][j]);
        nexStrm1 << std::endl;
        }
    nexStrm1 << "   ;" << std::endl;
    nexStrm1 << "end;" << std::endl << std::endl;
    
    nexStrm1 << "[MAXIMUMALLY PARTITIONED MIXED MODEL]" << std::endl;
    nexStrm1 << "begin mrbayes;" << std::endl;
	nexStrm1 << "   set autoclose=yes;" << std::endl;
	nexStrm1 << "   log start filename=sim_mb_max_" << repStr << "_log.txt;" << std::endl;
	begSite = 1;
    endSite = 0;
    for (int i=0; i<numSites.size(); i++)
        {
        endSite += numSites[i];
        nexStrm1 << "   charset part_" << i+1 << " = " << begSite << "-" << endSite << ";" << std::endl;
        begSite = endSite + 1;
        }
    nexStrm1 << "   partition my_part = " << numSites.size() << ":";
    for (int i=0; i<numSites.size(); i++)
        {
        nexStrm1 << "part_" << i+1;
        if ( i + 1 != numSites.size() )
            nexStrm1 << ",";
        }
    nexStrm1 << ";" << std::endl;
    nexStrm1 << "   set partition=my_part;" << std::endl;
    nexStrm1 << "   lset applyto=(all) nst=6 rates=gamma;" << std::endl;
    nexStrm1 << "   prset revmatpr=dirichlet(1,1,1,1,1,1) statefreqpr=dirichlet(1,1,1,1) shapepr=uniform(0.1,50);" << std::endl;
    nexStrm1 << "   unlink statefreq=(all) revmat=(all) shape=(all);" << std::endl;
	nexStrm1 << "   prset applyto=(all) ratepr=variable;" << std::endl;
	nexStrm1 << "   mcmc ngen=2000000 printfreq=1000 samplefreq=1000" << std::endl;
	nexStrm1 << "   nchains=1 savebrlens=yes filename=sim_mb_max_" << repStr << ";" << std::endl;
	nexStrm1 << "   sumt filename=sim_mb_max_" << repStr << " burnin=1000;" << std::endl;
	nexStrm1 << "   sump filename=sim_mb_max_" << repStr << " burnin=1000;" << std::endl;
	nexStrm1 << "end;" << std::endl;
	nexStrm1 << std::endl;
 
	// output the nexus file analysis under the true mixed model
    nexStrm2 << "#NEXUS" << std::endl << std::endl;
    nexStrm2 << "begin data;" << std::endl;
    nexStrm2 << "   dimensions ntax=" << settingsPtr->getNumTaxa() << " nchar=" << totalLength << ";" << std::endl;
    nexStrm2 << "   format datatype=dna;" << std::endl;
    nexStrm2 << "   matrix" << std::endl;
    for (int i=0; i<settingsPtr->getNumTaxa(); i++)
	{
        nexStrm2 << "   Taxon_" << i+1 << "   ";
        if (i+1 < 10)
            nexStrm2 << " ";
        for (int j=0; j<totalLength; j++)
            nexStrm2 << convertToNuc(matrix[i][j]);
        nexStrm2 << std::endl;
	}
    nexStrm2 << "   ;" << std::endl;
    nexStrm2 << "end;" << std::endl << std::endl;
    
    nexStrm2 << "[TRUE MIXED MODEL]" << std::endl;
	nexStrm2 << "begin mrbayes;" << std::endl;
	nexStrm2 << "   set autoclose=yes;" << std::endl;
	nexStrm2 << "   log start filename=sim_mb_tru_" << repStr << "_log.txt;" << std::endl;
	begSite = 1;
    endSite = 0;
    for (int i=0; i<numSites.size(); i++)
	{
        endSite += numSites[i];
        nexStrm2 << "   charset part_" << i+1 << " = " << begSite << "-" << endSite << ";" << std::endl;
        begSite = endSite + 1;
	}
    nexStrm2 << "   partition my_part = " << numSites.size() << ":";
    for (int i=0; i<numSites.size(); i++)
	{
        nexStrm2 << "part_" << i+1;
        if ( i + 1 != numSites.size() )
            nexStrm2 << ",";
	}
    nexStrm2 << ";" << std::endl;
    nexStrm2 << "   set partition=my_part;" << std::endl;
    nexStrm2 << "   lset applyto=(all) nst=6 rates=gamma;" << std::endl;
    nexStrm2 << "   prset revmatpr=dirichlet(1,1,1,1,1,1) statefreqpr=dirichlet(1,1,1,1) shapepr=uniform(0.1,50);" << std::endl;
    nexStrm2 << "   unlink revmat=(all) statefreq=(all) shape=(all);" << std::endl;
    nexStrm2 << "   link revmat=(1,2,3);" << std::endl;
    nexStrm2 << "   link revmat=(4,5,6);" << std::endl;
    nexStrm2 << "   link statefreq=(1,4);" << std::endl;
    nexStrm2 << "   link statefreq=(2,5);" << std::endl;
    nexStrm2 << "   link statefreq=(3,6);" << std::endl;
    nexStrm2 << "   link shape=(1,6);" << std::endl;
    nexStrm2 << "   link shape=(2,5);" << std::endl;
    nexStrm2 << "   link shape=(3,4);" << std::endl;
	nexStrm2 << "   prset applyto=(3,6) ratepr=variable;" << std::endl;
	nexStrm2 << "   mcmc ngen=2000000 printfreq=1000 samplefreq=1000" << std::endl;
	nexStrm2 << "   nchains=1 savebrlens=yes filename=sim_mb_tru_" << repStr << ";" << std::endl;
	nexStrm2 << "   sumt filename=sim_mb_tru_" << repStr << " burnin=1000;" << std::endl;
	nexStrm2 << "   sump filename=sim_mb_tru_" << repStr << " burnin=1000;" << std::endl;
	nexStrm2 << "end;" << std::endl;
	
	// close the file streams
    apStrm.close();
    nexStrm1.close();
    nexStrm2.close();
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