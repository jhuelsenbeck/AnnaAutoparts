#include <iostream>
#include <iomanip>
#include <set>
#include "Alignment.h"
#include "Chunk.h"
#include "MbMatrix.h"
#include "MbTransitionMatrix.h"
#include "Model.h"
#include "ParmAsrv.h"
#include "ParmFreqs.h"
#include "ParmLength.h"
#include "ParmSubrates.h"
#include "ParmTree.h"
#include "Settings.h"



Chunk::Chunk(Alignment* ap, Settings* sp, Model* mp, int si) {

	// remember the location of important objects
	alignmentPtr = ap;
	modelPtr     = mp;
	settingsPtr  = sp;
	
	// initialize variables
	subsetId     = si;
    update       = false;
	numNodes     = 2 * alignmentPtr->getNumTaxa() - 2;
	numGammaCats = settingsPtr->getNumGammaCats();
	
	// get the list of sites that are included in this "chunk", or subset
	for (int i=0; i<alignmentPtr->getNumChar(); i++)
		{
		if ( alignmentPtr->getIsExcluded(i) == false && alignmentPtr->getPartitionId(i) == subsetId + 1 )
			includedSites.insert(i);
		}
	numSites = includedSites.size();
	
	// allocate the space for the conditional likelihoods for the chunk
	int condLikeSize = numNodes * numSites * 4 * numGammaCats;
	cls = new double[condLikeSize];
	for (int i=0; i<condLikeSize; i++)
		cls[i] = 0.0;
	clsPtr = new double*[numNodes];
	for (int i=0; i<numNodes; i++)
		clsPtr[i] = &cls[i * numSites * 4 * numGammaCats];
		
	// initialize the conditional likelihoods for the chunk
	for (int i=0; i<alignmentPtr->getNumTaxa(); i++)
		{
		double* clp = clsPtr[i];
		for (std::set<int>::iterator c=includedSites.begin(); c != includedSites.end(); c++)
			{
			int nucCode = alignmentPtr->getNucleotide(i, (*c));
			int nucs[4];
			alignmentPtr->getPossibleNucs(nucCode, nucs);
			for (int k=0; k<numGammaCats; k++)
				{
				for (int s=0; s<4; s++)
					{
					if (nucs[s] == 1)
						clp[s] = 1.0;
					}
				clp += 4;
				}
			}
		}
		
	// allocate the transition probability matrices
	tis = new MbMatrix<double>*[numNodes];
	tis[0] = new MbMatrix<double>[numGammaCats*numNodes];
	for (int i=1; i<numNodes; i++)
		tis[i] = tis[i-1] + numGammaCats;
	for (int i=0; i<numNodes; i++)
		for (int k=0; k<numGammaCats; k++)
			tis[i][k] = MbMatrix<double>(4,4);

	// set up the transition probability calculator
	std::vector<double> bf(4, 0.25);
	std::vector<double> sr(6, 1.0);
	tiMatrix = new MbTransitionMatrix(sr, bf,  true);
		
#	if 0
	std::cout << "Sites for subset " << subsetId << ": ";
	for (std::set<int>::iterator p=includedSites.begin(); p != includedSites.end(); p++)
		std::cout << (*p) << " ";
	std::cout << std::endl;
	printTipCls();
#	endif
}

Chunk::~Chunk(void) {

	delete [] cls;
	delete [] clsPtr;
	delete [] tis[0];
	delete [] tis;
	delete tiMatrix;
}

double Chunk::lnLikelihood(bool storeScore) {

	// update the transition probabilities
	updateTransitionProbabilities();
	
	// get some variables that are used over-and-over again
	int stride = 4 * numGammaCats;
	Tree* treePtr = modelPtr->findTree(subsetId);
		
	/* pass down tree, filling in conditional likelihoods using the sum-product algorithm
	   (Felsenstein pruning algorithm) */
	double *lnScaler = new double[numSites];
	for (int i=0; i<numSites; i++)
		lnScaler[i] = 0.0;
		
	for (int n=0; n<treePtr->getNumNodes(); n++)
		{
		Node* p = treePtr->getDownPassNode( n );
		if ( p->getIsLeaf() == false )
			{
			if (p->getAnc()->getAnc() == NULL)
				{
				/* three-way split */
				int lftIdx = p->getLft()->getIndex();
				int rhtIdx = p->getRht()->getIndex();
				int ancIdx = p->getAnc()->getIndex();
				int idx    = p->getIndex();
				double *clL = getClsForNode(lftIdx);
				double *clR = getClsForNode(rhtIdx);
				double *clA = getClsForNode(ancIdx);
				double *clP = getClsForNode(idx   );
				for (int c=0; c<numSites; c++)
					{
					for (int k=0; k<numGammaCats; k++)
						{
						for (int i=0; i<4; i++)
							{
							double sumL = 0.0, sumR = 0.0, sumA = 0.0;
							for (int j=0; j<4; j++)
								{
								sumL += clL[j] * tis[lftIdx][k][i][j];
								sumR += clR[j] * tis[rhtIdx][k][i][j];
								sumA += clA[j] * tis[idx   ][k][i][j];
								}
							clP[i] = sumL * sumR * sumA;
							}
						clP += 4;
						clL += 4;
						clR += 4;
						clA += 4;
						}
					}
				}
			else
				{
				/* two-way split */
				int lftIdx = p->getLft()->getIndex();
				int rhtIdx = p->getRht()->getIndex();
				int idx    = p->getIndex();
				double *clL = getClsForNode(lftIdx);
				double *clR = getClsForNode(rhtIdx);
				double *clP = getClsForNode(idx   );
				for (int c=0; c<numSites; c++)
					{
					for (int k=0; k<numGammaCats; k++)
						{
						for (int i=0; i<4; i++)
							{
							double sumL = 0.0, sumR = 0.0;
							for (int j=0; j<4; j++)
								{
								sumL += clL[j] * tis[lftIdx][k][i][j];
								sumR += clR[j] * tis[rhtIdx][k][i][j];
								}
							clP[i] = sumL * sumR;
							}
						clP += 4;
						clL += 4;
						clR += 4;
						}
					}
				}
				
			/* scale */
#			if 1
			double *clP = getClsForNode( p->getIndex() );
			for (int c=0; c<numSites; c++)
				{
				double maxVal = 0.0;
				for (int i=0; i<stride; i++)
					{
					if (clP[i] > maxVal)
						maxVal = clP[i];
					}
				double scaler = 1.0 / maxVal;
				for (int i=0; i<stride; i++)
					clP[i] *= scaler;
				lnScaler[c] += log(maxVal);
				clP += stride;
				}
#			endif
						
			}
		}
		
	/* calculate likelihood */
	Node *p = treePtr->getRoot()->getLft();
	std::vector<double> f = modelPtr->findBaseFreqs(subsetId)->getFreq();
	double catProb = 1.0 / numGammaCats;
	double *clP = getClsForNode( p->getIndex() );
	double lnL = 0.0;
	for (int c=0; c<numSites; c++)
		{
		double like = 0.0;
		for (int k=0; k<numGammaCats; k++)
			{
			for (int i=0; i<4; i++)
				like += clP[i] * f[i] * catProb;
			clP += 4;
			}
		lnL += log( like ) + lnScaler[c];
		}
		
	delete [] lnScaler;
	
    //std::cout << "lnL[" << subsetId << "] = " << lnL << std::endl;
    
    if (storeScore == true)
        storedLnL = lnL;
    update = false;
    mostRecentLnL = lnL;
    
	return lnL;
}

void Chunk::printTis(void) {

	for (int nde=0; nde<numNodes; nde++)
        {
        std::cout << "Transition probabilities for node " << nde << std::endl;
        for (int i=0; i<4; i++)
            {
            std::cout << "   ";
            for (int k=0; k<numGammaCats; k++)
                {
                for (int j=0; j<4; j++)
                    std::cout << std::fixed << std::setprecision(4) << tis[nde][k][i][j] << " ";
                if (k + 1 != numGammaCats)
                    std::cout << " -- ";
                }
            std::cout << std::endl;
            }
        }
}

void Chunk::printTipCls(void) {

	for (int c=0; c<numSites; c++)
		{
		for (int k=0; k<numGammaCats; k++)
			{
			if (k == 0)
				std::cout << std::setw(4) << c << " -- ";
			else 
				std::cout << "        ";
			for (int i=0; i<alignmentPtr->getNumTaxa(); i++)
				{
				double *clp = clsPtr[i] + (c * 4 * numGammaCats) + 4*k;
				for (int s=0; s<4; s++)
					std::cout << std::fixed << std::setprecision(0) << clp[s];
				std::cout << " ";
				}
			std::cout << std::endl;
			}
		}
}

void Chunk::updateTransitionProbabilities(void) {

	// get poiters to the parameters that are necessary to properly update the transition probabilities
	Tree*       treePtr       = modelPtr->findTree(subsetId);
	TreeLength* treeLengthPtr = modelPtr->findTreeLength(subsetId);
	SubRates*   subRatesPtr   = modelPtr->findSubRates(subsetId);
	BaseFreqs*  baseFreqsPtr  = modelPtr->findBaseFreqs(subsetId);
	Asrv*       asrvPtr       = modelPtr->findAsrv(subsetId);
	
	// get the parameter values
	std::vector<double> bf = baseFreqsPtr->getFreq();
	std::vector<double> sr = subRatesPtr->getSubRate();
	std::vector<double> r  = asrvPtr->getRate();
	double len             = treeLengthPtr->getLength();
    
	// update the rate matrix (also updates the eigensystem
	tiMatrix->updateQ(sr, bf);

	// update the transition probabilities
	for (int n=0; n<treePtr->getNumNodes(); n++)
		{
		Node *p = treePtr->getDownPassNode( n );
		if ( p->getAnc() != NULL )
			{
			int idx = p->getIndex();
			double prop = p->getP();
			double v    = len * prop;
			for (int k=0; k<numGammaCats; k++)
				{
				double vr = v * r[k];
				tis[idx][k] = tiMatrix->tiProbs(vr, tis[idx][k]);
				//std::cout << k << " " << vr << std::endl;
				//cout << ti[idx][k] << endl;
				}
			}
		}
        
#   if 0
    std::cout << "bf: ";
    for (int i=0; i<bf.size(); i++)
        std::cout << bf[i] << " ";
    std::cout << std::endl;
    std::cout << "sr: ";
    for (int i=0; i<sr.size(); i++)
        std::cout << sr[i] << " ";
    std::cout << std::endl;
    std::cout << "r: ";
    for (int i=0; i<r.size(); i++)
        std::cout << r[i] << " ";
    std::cout << std::endl;
    std::cout << "len: " << len << std::endl;
#   endif
}



