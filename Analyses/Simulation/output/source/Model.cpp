#include "Alignment.h"
#include "Chunk.h"
#include "MbRandom.h"
#include "Model.h"
#include "ParmAsrv.h"
#include "ParmFreqs.h"
#include "ParmLength.h"
#include "ParmSubrates.h"
#include "ParmTree.h"
#include "Restaurant.h"
#include "StateSets.h"
#include "Table.h"

#undef DEBUG_LNL
#define SMART_LIKES



Model::Model(Settings* sp, MbRandom* rp, Alignment* ap) {

	// remember the location of important objects
	settingsPtr  = sp;
	ranPtr       = rp;
	alignmentPtr = ap;
	
	// initialize important variables
	numSubsets   = alignmentPtr->getNumSubsets();

	// set up "chunks", containing the data for each subset
	for (int i=0; i<numSubsets; i++)
		chunks.push_back( new Chunk(alignmentPtr, settingsPtr, this, i) );
		
	// initialize and set up the state sets for calculating parsimony scores
	stateSets = new StateSets(alignmentPtr, settingsPtr);

	// set up the parameters
	restaurants.push_back( new Restaurant(ranPtr, settingsPtr, alignmentPtr, this, TREE,     numSubsets, 0.0) );
	restaurants.push_back( new Restaurant(ranPtr, settingsPtr, alignmentPtr, this, ASRV,     numSubsets, 0.5) );
	restaurants.push_back( new Restaurant(ranPtr, settingsPtr, alignmentPtr, this, SUBRATE,  numSubsets, 0.5) );
	restaurants.push_back( new Restaurant(ranPtr, settingsPtr, alignmentPtr, this, BASEFREQ, numSubsets, 0.5) );
	restaurants.push_back( new Restaurant(ranPtr, settingsPtr, alignmentPtr, this, LENGTH,   numSubsets, 0.5) );

	// calculate the probability of changing each parameter (note that the probabilities are in the same order as the restaurants, above)
	proposalProb.push_back( 3 );
	proposalProb.push_back( 1 );
	proposalProb.push_back( 1 );
	proposalProb.push_back( 1 );
	proposalProb.push_back( 1 );
	double sum = 0.0;
	for (int i=0; i<proposalProb.size(); i++)
		sum += proposalProb[i];
	for (int i=0; i<proposalProb.size(); i++)
		proposalProb[i] /= sum;
}

Model::~Model(void) {

	for (int i=0; i<numSubsets; i++)
		{
		delete chunks[i];
		}
	delete stateSets;
}

Asrv* Model::findAsrv(int part) {

	for (std::vector<Restaurant*>::iterator p=restaurants.begin(); p != restaurants.end(); p++)
		{
		if ( (*p)->getParmId() == ASRV )
			{
			Table* tbl = (*p)->getTableWithPatron(part);
			Asrv* theDerivedPtr = dynamic_cast<Asrv *>(tbl->getParm());
			if ( theDerivedPtr == 0  )
				{
				std::cerr << "ERROR: Problem downcasting ASRV parameter" << std::endl;
				exit(1);
				}
			return theDerivedPtr;
			}
		}
	std::cerr << "ERROR: Could not find ASRV parameter" << std::endl;
	return NULL;
	
}

BaseFreqs* Model::findBaseFreqs(int part) {

	for (std::vector<Restaurant *>::iterator p=restaurants.begin(); p != restaurants.end(); p++)
		{
		if ( (*p)->getParmId() == BASEFREQ )
			{
			Table* tbl = (*p)->getTableWithPatron(part);
			BaseFreqs* theDerivedPtr = dynamic_cast<BaseFreqs *>(tbl->getParm());
			if ( theDerivedPtr == 0  )
				{
				std::cerr << "ERROR: Problem downcasting BaseFreqs parameter" << std::endl;
				exit(1);
				}
			return theDerivedPtr;
			}
		}
	std::cerr << "ERROR: Could not find BaseFreqs parameter" << std::endl;
	return NULL;
	
}

SubRates* Model::findSubRates(int part) {

	for (std::vector<Restaurant *>::iterator p=restaurants.begin(); p != restaurants.end(); p++)
		{
		if ( (*p)->getParmId() == SUBRATE )
			{
			Table* tbl = (*p)->getTableWithPatron(part);
			SubRates* theDerivedPtr = dynamic_cast<SubRates *>(tbl->getParm());
			if ( theDerivedPtr == 0  )
				{
				std::cerr << "ERROR: Problem downcasting SubRates parameter" << std::endl;
				exit(1);
				}
			return theDerivedPtr;
			}
		}
	std::cerr << "ERROR: Could not find SubRates parameter" << std::endl;
	return NULL;
	
}

Tree* Model::findTree(int part) {

	for (std::vector<Restaurant *>::iterator p=restaurants.begin(); p != restaurants.end(); p++)
		{
		if ( (*p)->getParmId() == TREE )
			{
			Table* tbl = (*p)->getTableWithPatron(part);
			Tree* theDerivedPtr = dynamic_cast<Tree *>(tbl->getParm());
			if ( theDerivedPtr == 0  )
				{
				std::cerr << "ERROR: Problem downcasting Tree parameter" << std::endl;
				exit(1);
				}
			return theDerivedPtr;
			}
		}
	std::cerr << "ERROR: Could not find Tree parameter" << std::endl;
	return NULL;
	
}

TreeLength* Model::findTreeLength(int part) {

	for (std::vector<Restaurant *>::iterator p=restaurants.begin(); p != restaurants.end(); p++)
		{
		if ( (*p)->getParmId() == LENGTH )
			{
			Table* tbl = (*p)->getTableWithPatron(part);
			TreeLength* theDerivedPtr = dynamic_cast<TreeLength *>(tbl->getParm());
			if ( theDerivedPtr == 0  )
				{
				std::cerr << "ERROR: Problem downcasting TreeLength parameter" << std::endl;
				exit(1);
				}
			return theDerivedPtr;
			}
		}
	std::cerr << "ERROR: Could not find TreeLength parameter" << std::endl;
	return NULL;
}

int Model::getRestaurantId(Restaurant* rest) {

	for (int i=0; i<restaurants.size(); i++)
		{
		if (restaurants[i] == rest)
			return i;
		}
	return -1;
}

Restaurant* Model::getRestaurantToChange(void) {

	double u = ranPtr->uniformRv();
	double sum = 0.0;
	int whichRestaurant = 0;
	for (int i=0; i<proposalProb.size(); i++)
		{
		sum += proposalProb[i];
		if (u < sum)
			{
			whichRestaurant = i;
			break;
			}
		}
	return restaurants[whichRestaurant];
}

double Model::lnLikelihood(void) {

#	if defined(DEBUG_LNL)
	return 0.0;
#	endif
	double lnL = 0.0;
	for (int i=0; i<numSubsets; i++)
		{
#		if defined (SMART_LIKES)
        if ( chunks[i]->getUpdate() == true )
            lnL += chunks[i]->lnLikelihood(false);
        else 
            lnL += chunks[i]->getStoredLnL();
#		else
		lnL += chunks[i]->lnLikelihood(false);
#		endif
		}
	return lnL;
}

double Model::lnLikelihood(bool storeScore) {

#	if defined(DEBUG_LNL)
	return 0.0;
#	endif
	double lnL = 0.0;
	for (int i=0; i<numSubsets; i++)
		{
#		if defined (SMART_LIKES)
		lnL += chunks[i]->lnLikelihood(storeScore);
#		else
		lnL += chunks[i]->lnLikelihood(false);
#		endif
		}
    //std::cout << "lnL[] = " << lnL << std::endl;
	return lnL;
}

double Model::lnLikelihood(int patron, bool storeScore) {

#	if defined(DEBUG_LNL)
	return 0.0;
#	endif
#	if defined (SMART_LIKES)
	return chunks[patron]->lnLikelihood(storeScore);
#	else
	return chunks[patron]->lnLikelihood(false);
#	endif
}

void Model::setStoredLnLToMostRecentLnL(Table* tbl) {

   std::set<int> tablePatrons = tbl->getPatrons();
   for (std::set<int>::iterator p = tablePatrons.begin(); p != tablePatrons.end(); p++)
        chunks[(*p)]->setStoredLnLToMostRecentLnL();
}

void Model::setUpdateFlags(Table* tbl) {

    for (int i=0; i<numSubsets; i++)
        chunks[i]->setUpdate(false);
   std::set<int> tablePatrons = tbl->getPatrons();
   for (std::set<int>::iterator p = tablePatrons.begin(); p != tablePatrons.end(); p++)
        chunks[(*p)]->setUpdate(true);
}



