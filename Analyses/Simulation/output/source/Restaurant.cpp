#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include "Chunk.h"
#include "MbRandom.h"
#include "Model.h"
#include "Parm.h"
#include "Partition.h"
#include "Restaurant.h"
#include "Settings.h"
#include "Table.h"

#undef DEBUG_NORMALIZE



Restaurant::Restaurant(MbRandom* rp, Settings* sp, Alignment* ap, Model* mp, ParmId pid, int np, double ps) {

	// remember the location of some objects
	ranPtr            = rp;
	settingsPtr       = sp;
	alignmentPtr      = ap;
	modelPtr          = mp;
	parmId            = pid;

	// initalize parameters
	numPatrons        = np;           
	numAuxTables      = 5;  
	probUpdateSeating = ps; 

	// initialize the concentration parameter
	if (settingsPtr->getIsConcFixed() == true)
		{
		alpha = calcAlphaFromEt( settingsPtr->getExpNumCats() );
		std::cout << "alpha fixed to " << alpha << std::endl;
		}
	else
		{
		double m = settingsPtr->getPriorConcMean();
		double v = settingsPtr->getPriorConcVariance();
		gammaBeta  = m / v;
		gammaAlpha = gammaBeta * m;
		alpha = ranPtr->gammaRv(gammaAlpha, gammaBeta);
		}
		
	// set the name of the restaurant
	if (parmId == ASRV)
		{
		name = "asrv";
		shortName = "rv";
		}
	else if (parmId == TREE)
		{
		name = "tree";
		shortName = "t";
		}
	else if (parmId == SUBRATE)
		{
		name = "subrates";
		shortName = "sr";
		}
	else if (parmId == BASEFREQ)
		{
		name = "basefreq";
		shortName = "bf";
		}
	else if (parmId == LENGTH)
		{
		name = "treelength";
		shortName = "tl";
		}
		
	// seat the patrons in the restaurant
	if (probUpdateSeating < 0.0000001)
		{
		Table* tbl = new Table( ranPtr, settingsPtr, alignmentPtr, modelPtr, parmId );
		tables.insert( tbl );
		for (int i=0; i<numPatrons; i++)
			tbl->seatPatron( i );
		}
	else 
		{
#		if 1
		for (int i=0; i<numPatrons; i++)
			{
			double newTableProb = alpha / (i + alpha);
			double u = ranPtr->uniformRv();
			if (u < newTableProb)
				{
				Table* tbl = new Table( ranPtr, settingsPtr, alignmentPtr, modelPtr, parmId );
				tables.insert( tbl );
				tbl->seatPatron( i );
				}
			else 
				{
				Table* tbl = pickTableAtRandomFromPrior();
				tbl->seatPatron( i );
				}
			}
#		else
		for (int i=0; i<numPatrons; i++)
			{
			Table* tbl = new Table( ranPtr, settingsPtr, alignmentPtr, modelPtr, parmId );
			tables.insert( tbl );
			tbl->seatPatron( i );
			}
		probUpdateSeating = 0.0;
#		endif
		}

	// open files for output
	std::string parmFileName = settingsPtr->getOutPutFileName() + "." + name + ".out";
	std::string partFileName = settingsPtr->getOutPutFileName() + "." + name + ".part";
	parmOut.open( parmFileName.c_str(), std::ios::out );
	partOut.open( partFileName.c_str(), std::ios::out );
	if ( probUpdateSeating > 0.0 )
		{
		parmOut << "Generation" << '\t';
		parmOut << "ConcentrationParm" << '\t';
		parmOut << "Degree" << '\t';
		partOut << "Degree" << '\t';
		partOut << "Partition" << '\t';
		}
	partOut << std::endl;
	if (probUpdateSeating == 0.0)
		{
		std::string headerStr = getTableWithPatron(0)->getParm()->getParmHeader(-1);
		parmOut << headerStr;
		}
	else
		{
		for (int i=0; i<numPatrons; i++)
			{
			std::string headerStr = getTableWithPatron(0)->getParm()->getParmHeader(i+1);
			parmOut << headerStr;
			}
		}
	parmOut << std::endl;
	
	// allocate the RGF
	rgf = new int[numPatrons];
	
	print();
}

Restaurant::~Restaurant(void) {

	for (std::set<Table *>::iterator p=tables.begin(); p != tables.end(); p++)
		delete (*p);
	delete [] rgf;
	parmOut.close();
}

void Restaurant::addPartition(void) {

	setRgf();
	Partition* aPartition = new Partition(numPatrons, rgf);
	sampledPartitions.push_back( aPartition );
}

double Restaurant::acceptanceProb(double lnR) {

	if (lnR < -300.0)
		return 0.0;
	else if (lnR > 0.0)
		return 1.0;
	else
		return exp( lnR );
}

double Restaurant::calcAlphaFromEt(double expT) {

	double a = 0.000001;
	double ea = expNumTables(a);
	bool goUp;
	if (ea < expT)
		goUp = true;
	else
		goUp = false;
	double increment = 0.1;
	while ( fabs(ea - expT) > 0.000001 )
		{
		if (ea < expT && goUp == true)
			{
			a += increment;
			}
		else if (ea > expT && goUp == false)
			{
			a -= increment;
			}
		else if (ea < expT && goUp == false)
			{
			increment /= 2.0;
			goUp = true;
			a += increment;
			}
		else
			{
			increment /= 2.0;
			goUp = false;
			a -= increment;
			}
		ea = expNumTables(a);
		}
	//cout << ea << " <-> " << expT << " " << "alpha=" << a << endl;
	return a;
}

bool Restaurant::change(void) {

	bool wasAccepted = false;
	double u = ranPtr->uniformRv();
	if (u < probUpdateSeating)
		{
		for (int i=0; i<numPatrons; i++)
			updateSeating(i);
		wasAccepted = true;
		for (std::set<Table*>::iterator tbl = tables.begin(); tbl != tables.end(); tbl++)
			changeParmOnTable(*tbl);
		}
	else
		{
		Table* tbl = pickTableUniformlyAtRandom();
		wasAccepted = changeParmOnTable( tbl );
		}

	// update the concentration parameter
	if ( settingsPtr->getIsConcFixed() == false && probUpdateSeating > 0.000001 )
		alpha = sampleAlpha( getNumTables(), numPatrons, alpha, gammaAlpha, gammaBeta );
	
	return wasAccepted;
}

bool Restaurant::changeParmOnTable(Table* tbl) {

	Parm* oldParm = tbl->getParm();                     // get a pointer to the other (copy) parameter
	tbl->flipCurParmId();                               // flip the parm id so we can...
	Parm* newParm = tbl->getParm();                     // get a pointer to the parameter on the table
	double oldLnPrior = newParm->lnPriorProb();         // get the prior probability of the old state of the parameter
    double oldLnL = modelPtr->lnLikelihood();
	double lnProposalProb = newParm->update();          // update the parameter, capturing the log of the Hastings ratio
	double newLnPrior = newParm->lnPriorProb();         // calculate the prior probability of the updated state of the parameter...
	double lnPriorRatio = newLnPrior - oldLnPrior;      // so we can calculate the log of the prior ratio
    
    modelPtr->setUpdateFlags(tbl);
    double newLnL = modelPtr->lnLikelihood();
	double lnLikelihoodRatio = newLnL - oldLnL;         // allow us to calculate the log of the likelihood ratio
	double r = acceptanceProb( lnLikelihoodRatio + lnPriorRatio + lnProposalProb );
	
	bool acceptUpdate = false;
	if ( ranPtr->uniformRv() < r )
		acceptUpdate = true;
	if (acceptUpdate == true)
		{
        modelPtr->setStoredLnLToMostRecentLnL(tbl);
		(*oldParm) = (*newParm);
		}
	else
		{
		(*newParm) = (*oldParm);
		}
	
	return acceptUpdate;
}

void Restaurant::deleteUnoccupiedTables(void) {

	std::vector<Table*> tablesToDelete;
	for (std::set<Table*>::iterator p = tables.begin(); p != tables.end(); p++)
		{
		if ( (*p)->numPatronsAtTable() == 0 )
			tablesToDelete.push_back( (*p) );
		}
	for (std::vector<Table*>::iterator p = tablesToDelete.begin(); p != tablesToDelete.end(); p++)
		removeTable(*p);
}

double Restaurant::expNumTables(double a) {

	double expectedNum = 0.0;
	for (int i=1; i<=numPatrons; i++)
		expectedNum += ( 1.0 / (i - 1.0 + a) );
	expectedNum *= a;
	return expectedNum;
}

Table* Restaurant::getTableWithPatron(int patron) {

	for (std::set<Table*>::iterator p = tables.begin(); p != tables.end(); p++)
		{
		if ( (*p)->isPatronAtTable(patron) == true )
			return (*p);
		}
	return NULL;
}

void Restaurant::normalizeLogProbabilitiesInMap(std::map<Table*,double>& m) {

#	if defined(DEBUG_NORMALIZE)
	int i = 0;
	for (std::map<Table*,double>::iterator p = m.begin(); p != m.end(); p++)
		std::cout << i++ << " -- " << p->second << std::endl;
#	endif

	// set up an iterator for the map
	std::map<Table*,double>::iterator p = m.begin();
	
	// find the largest value in the set of log probabilities
	double lnC = p->second;
	for (p = m.begin(); p != m.end(); p++)
		{
		if (p->second > lnC)
			lnC = p->second;
		}
	
	// pull the log factor lnC out of every log probability in the set
	for (p = m.begin(); p != m.end(); p++)
		p->second -= lnC;

	// exponentiate the log probabilities to make them probabilities
	double sum = 0.0;
	for (p = m.begin(); p != m.end(); p++)
		{
		if ( p->second < -300.0 )
			p->second = 0.0;
		else
			p->second = exp( p->second );
		sum += p->second;
		}
	
	// normalize
	for (p = m.begin(); p != m.end(); p++)
		p->second /= sum;
		
#	if defined(DEBUG_NORMALIZE)
	i = 0;
	for (p = m.begin(); p != m.end(); p++)
		std::cout << i++ << " -- " << p->second << std::endl;
#	endif

}

Table* Restaurant::pickTableAtRandomFromPrior(void) {

	// get the number of patrons currently seated at all of the tables
	int n = 0;
	for (std::set<Table*>::iterator p = tables.begin(); p != tables.end(); p++)
		n += (*p)->numPatronsAtTable();
		
	double u = ranPtr->uniformRv();
	double sum = 0.0;
	for (std::set<Table*>::iterator p = tables.begin(); p != tables.end(); p++)
		{
		sum += (double)((*p)->numPatronsAtTable()) / n;
		if (u < sum)
			return (*p);
		}
	return NULL;
}

Table* Restaurant::pickTableUniformlyAtRandom(void) {

	int whichTable = (int)(ranPtr->uniformRv() * tables.size());
	int i = 0;
	for (std::set<Table*>::iterator p = tables.begin(); p != tables.end(); p++)
		{
		if (i == whichTable)
			return (*p);
		i++;
		}
	return NULL;
}

void Restaurant::print(void) {

	int i = 0;
	for (std::set<Table*>::iterator p=tables.begin(); p != tables.end(); p++)
		{
		std::cout << "Table " << ++i << std::endl;
		(*p)->print();
		}
}

void Restaurant::removeTable(Table* tbl) {

	tables.erase( tbl );
	delete tbl;
}

double Restaurant::sampleAlpha(int k, int n, double oldAlpha, double a, double b) {

	/* Step 1: Draw a Beta(oldAlpha+1, n) distribution to get eta */
	std::vector<double> z(2);
	std::vector<double> f(2);
	z[0] = oldAlpha + 1.0;
	z[1] = (double)n;
	ranPtr->dirichletRv(z, f);
	double eta = f[0];
	
	/* Step 2: Draw a new value for alpha, based on k and eta */
	double u = ranPtr->uniformRv();
	double x = ( a + (double)k - 1.0 ) / ( (double)n * (b - log(eta)) );
	double newAlpha;
	if ( (u / (1.0 - u)) < x)
		newAlpha = ranPtr->gammaRv(a + k, b - log(eta));
	else
		newAlpha = ranPtr->gammaRv(a + k - 1.0, b - log(eta));
		
	return newAlpha;
}

void Restaurant::saveState(int n) {

	if ( probUpdateSeating > 0.0 )
		{
		parmOut << n << '\t';
		parmOut << std::fixed << std::setprecision(6) << alpha << '\t';
		parmOut << getNumTables() << '\t';
		}
	partOut << getNumTables() << '\t';
	setRgf();
	partOut << "[";
	for (int i=0; i<numPatrons; i++)
		{
		partOut << rgf[i];
		if (i + 1 < numPatrons)
			partOut << ",";
		}
	partOut << "] \t";
		
	if ( probUpdateSeating == 0.0 )
		{
		Table* tbl = getTableWithPatron( 0 );
		std::string temp = tbl->getParmString(n);
		parmOut << temp;
		}
	else
		{
		for (int i=0; i<numPatrons; i++)
			{
			Table* tbl = getTableWithPatron( i );
			std::string temp = tbl->getParmString(n);
			parmOut << temp;
			}
		}
		
	parmOut << std::endl;
	partOut << std::endl;
}

void Restaurant::setRgf(void) {

	for (int i=0; i<numPatrons; i++)
		rgf[i] = -1;
	int k = 1;
	for (int i=0; i<numPatrons; i++)
		{
		if (rgf[i] == -1)
			{
			Table* tbl = getTableWithPatron( i );
			std::set<int> tablePatrons = tbl->getPatrons();
			for (std::set<int>::iterator elem = tablePatrons.begin(); elem != tablePatrons.end(); elem++)
				{
				rgf[*elem] = k;
				}
			k++;
			}
		}
}

void Restaurant::simulateFromPrior(std::vector<int>& x, int nReps) {

	x.resize(numPatrons);
	for (int i=0; i<x.size(); i++)
		x[i] = 0;
	for (int i=0; i<nReps; i++)
		{
		double a = alpha;
		if (settingsPtr->getIsConcFixed() == false)
			a = ranPtr->gammaRv(gammaAlpha, gammaBeta);
		int nt = 0;
		for (int j=0; j<numPatrons; j++)
			{
			double pNew = a / (j + a);
			double u = ranPtr->uniformRv();
			if (u < pNew)
				nt++;
			}
		x[nt-1]++;
		}
#	if 0
	for (int i=0; i<numPatrons; i++)
		std::cout << std::setw(4) << i << "  " << x[i] << std::endl;
#	endif
}

void Restaurant::summarizePartitions(void) {

	// print the restaurant name
	std::cout << "Restaurant: " << getName() << std::endl;
	// first, calculate the mean partition
	Partition avePart(sampledPartitions);
	std::cout << "   Average Partition = ";
	avePart.print();
	
	// get the list of unique partitions
	vector<int> numOfThisPart;
	vector<Partition*> uniqueParts;
	vector<int> indices;
	int j = 0;
	for (vector<Partition*>::iterator p=sampledPartitions.begin(); p != sampledPartitions.end(); p++)
		{
		bool foundPart = false;
		int i = 0;
		for (vector<Partition*>::iterator q=uniqueParts.begin(); q != uniqueParts.end(); q++)
			{
			if ( (**p) == (**q) )
				{
				foundPart = true;
				break;
				}
			i++;
			}
		if ( foundPart == false )
			{
			uniqueParts.push_back( new Partition(**p) );
			numOfThisPart.push_back( 1 );
			indices.push_back( j++ );
			}
		else
			{
			numOfThisPart[i]++;
			}
		}
	mySort( numOfThisPart, indices, numOfThisPart.size() );
		
	std::cout << "   Unique partitions:" << endl;
	double cumulative = 0.0;
	for (int i=0; i<uniqueParts.size(); i++)
		{
		double prob = (double)numOfThisPart[i] / sampledPartitions.size();
		cumulative += prob;
		std::cout << "   " << std::setw(5) << i + 1 << " -- " << setw(5) << numOfThisPart[i] << " " << fixed << setprecision(4) << prob << " " << cumulative << " ";
		uniqueParts[ indices[i] ]->print();
		}
	
	// calculate the probability distribution for the number of sampled tables
	std::vector<int> px;
	int nReps = 10000;
	simulateFromPrior(px, nReps);

	std::cout << "   Posterior and prior probability distributions for the number of tables:" << std::endl;
	std::vector<int> pp(numPatrons+1, 0);
	int largestDegree = 0, n = 0;
	for (std::vector<Partition *>::iterator p=sampledPartitions.begin(); p != sampledPartitions.end(); p++)
		{
		pp[ (*p)->getDegree() ]++;
		if ( (*p)->getDegree() > largestDegree )
			largestDegree = (*p)->getDegree();
		n++;
		}
	for (int i=1; i<=numPatrons; i++)
		std::cout << "   " << std::setw(5) << i << " -- " << std::setw(5) << pp[i] << " " << std::fixed << std::setprecision(4) << (double)pp[i] / n << " " << std::fixed << std::setprecision(4) << (double)px[i-1] / nReps << std::endl;
	
	// calculate the probability that patrons are grouped together at the same table
	int** togetherness = new int*[numPatrons];
	togetherness[0] = new int[numPatrons*numPatrons];
	for (int i=1; i<numPatrons; i++)
		togetherness[i] = togetherness[i-1] + numPatrons;
	for (int i=0; i<numPatrons; i++)
		for (int j=0; j<numPatrons; j++)
			togetherness[i][j] = 0;
			
	for (std::vector<Partition *>::iterator p=sampledPartitions.begin(); p != sampledPartitions.end(); p++)
		{
		for (int i=0; i<numPatrons; i++)
			{
			for (int j=i+1; j<numPatrons; j++)
				{
				if ( (*p)->getElement(i) == (*p)->getElement(j) )
					togetherness[i][j]++;
				}
			}
		}
	
	std::cout << "   Probability that patrons are seated at the same table:" << std::endl;
	for (int i=0; i<numPatrons; i++)
		{
		for (int j=i+1; j<numPatrons; j++)
			{
			std::cout << "   " << std::setw(4) << i+1 << " " << std::setw(4) << j+1 << " -- " << std::fixed << std::setprecision(4) << (double)togetherness[i][j] / sampledPartitions.size() << std::endl;
			}
		}
		
	delete [] togetherness[0];
	delete [] togetherness;
}

void Restaurant::mySort(std::vector<int> &item, std::vector<int> &assoc, int count) {

	sort2(item, assoc, 0, count-1);
}

void Restaurant::sort2(std::vector<int> &item, std::vector<int> &assoc, int left, int right) {

	int i = left;
	int j = right;
	int x = item[(left+right)/2];
	do 
		{
		/*while (item[i] < x && i < right)
			i++;
		while (x < item[j] && j > left)
			j--;*/

		while (item[i] > x && i < right)
			i++;
		while (x > item[j] && j > left)
			j--;

		if (i <= j)
			{
			int y = item[i];
			item[i] = item[j];
			item[j] = y;
			
			int temp = assoc[i];
			assoc[i] = assoc[j];
			assoc[j] = temp;
			
			i++;
			j--;
			}
		} while (i <= j);
	if (left < j)
		sort2 (item, assoc, left, j);
	if (i < right)
		sort2 (item, assoc, i, right);

}

void Restaurant::updateSeating(int patron) {

	// remove patron from its current table 
	Table* tbl = getTableWithPatron( patron );
	tbl->removePatron( patron );
	if ( tbl->numPatronsAtTable() <= 0 )
		removeTable( tbl );
		
	// make some auxiliary tables
	for (int i=0; i<numAuxTables; i++)
		{
		Table* tbl = new Table(ranPtr, settingsPtr, alignmentPtr, modelPtr, parmId);
		tables.insert( tbl ); 
		}

	// make a map that contains as the key a pointer to the table and as the value the probability of seating the patron at that table
	std::map<Table*,double> probs, lnLikes;
	for (std::set<Table*>::iterator p = tables.begin(); p != tables.end(); p++)
		probs.insert( std::make_pair((*p),0.0) );
		
	// calculate the log probability (up to a constant) of seating the patron at each of the possible tables 
	for (std::map<Table*,double>::iterator p = probs.begin(); p != probs.end(); p++)
		{
		Table* tbl = p->first;                                    // get a pointer to the table
		tbl->seatPatron(patron);                                  // seat the patron at the i-th existing table
		double lnLPart = modelPtr->lnLikelihood(patron, false);   // calculate the likelihood when the patron is seated at the i-th table
		//double lnLPart = modelPtr->lnLikelihood(false);   // calculate the likelihood when the patron is seated at the i-th table
        p->second = lnLPart;
        lnLikes.insert( std::make_pair(tbl,lnLPart) );
		if (tbl->numPatronsAtTable() == 1)
			p->second += log(alpha / numAuxTables);
		else 
			p->second += log( tbl->numPatronsAtTable() - 1 );
		tbl->removePatron( patron );
		}


	// calculate the probabilities of reseating the element at each of the tables, and pick a table to reseat the element at
	normalizeLogProbabilitiesInMap(probs);
	
	// reseat the patron at one of the tables chosen in proportion to the probabilities calculated above
	double u = ranPtr->uniformRv();
	double sum = 0.0;
    Table* newTable = NULL;
	for (std::map<Table*,double>::iterator p = probs.begin(); p != probs.end(); p++)
		{
		sum += p->second;
		if (u < sum)
			{
			p->first->seatPatron(patron);
            newTable = p->first;
			break;
			}
		}
	
	// update the likelihood for that element
#	if 0
    std::map<Table*,double>::iterator it = lnLikes.find(newTable);
    if (it != lnLikes.end())
        modelPtr->getChunk(patron)->setStoredLnL(it->second);
    else 
        std::cerr << "Error: Problem storing log likelihood" << std::endl;
#	else 
	double lnLPart = modelPtr->lnLikelihood(patron, false);  
	modelPtr->getChunk(patron)->setStoredLnL(lnLPart);
#	endif

	// delete unoccupied tables 
	deleteUnoccupiedTables();
}


