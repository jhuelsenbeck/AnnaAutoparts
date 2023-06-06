#include <cmath>
#include <iomanip>
#include <iostream>
#include "Chunk.hpp"
#include "Model.hpp"
#include "Msg.hpp"
#include "Parameter.hpp"
#include "ParameterBaseFrequencies.hpp"
#include "ParameterExchangabilityRates.hpp"
#include "ParameterGammaShape.hpp"
#include "ParameterTree.hpp"
#include "ParameterTreeLength.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "Restaurant.hpp"
#include "Table.hpp"
#include "TableFactory.hpp"
#include "UpdateInfo.hpp"
#include "UserSettings.hpp"



Restaurant::Restaurant(const Restaurant& r) {

    modelPtr = r.modelPtr;
    settingsPtr = r.settingsPtr;
    alpha = r.alpha;
    numPatrons = r.numPatrons;
    isSeatingRv = r.isSeatingRv;
    gammaAlpha = r.gammaAlpha;
    gammaBeta = r.gammaBeta;
        
    parameter = copyParameter(r.parameter);
    if (parameter == NULL)
        Msg::error("Could not find parameter to copy in restaurant copy constructor");
        
    TableFactory& tf = TableFactory::tableFactory();
    for (Table* t : r.tables)
        {
        Table* newT = tf.getTable();
        newT->setRestaurant(this);
        std::set<int>& tPatrons = t->getPatrons();
        for (int i : tPatrons)
            newT->addPatron(i);
        newT->addParameter( copyParameter(t->getParameter()) );
        tables.insert(newT);
        }
}

Restaurant::Restaurant(Model* mp, UserSettings* s, bool sf, double a, int np, Parameter* parm) {

    isSeatingRv = sf;
    modelPtr    = mp;
    settingsPtr = s;
    numPatrons  = np;
    parameter   = parm;
    
    if (isSeatingRv == false)
        {
        // seating is not a random variable, so we seat all patrons at one table
        std::cout << "   * Parameter " << parameter->getName() << " subset arrangement is fixed (all subsets grouped together)" << std::endl;
        alpha = 0.0;
        Table* t = addTable();
        for (int i=0; i<numPatrons; i++)
            t->addPatron(i);
        std::cout << "   *    Fixed RGF = " << rgf() << std::endl;
        }
    else
        {
        // seating is a random variable following a DPP
        RandomVariable& rng = RandomVariable::randomVariableInstance();
        
        if (settingsPtr->getIsConcentrationParameterFixed() == true)
            {
            std::cout << "   * Parameter " << parameter->getName() << " subset arrangement is a random variable following a DPP with fixed alpha" << std::endl;
            alpha = a; // a is passed into the function based on E(#Tables)
            std::cout << "   *    Fixed Alpha = " << alpha << std::endl;
            }
        else
            {
            std::cout << "   * Parameter " << parameter->getName() << " subset arrangement is a random variable following a DPP with a gamma hyperprior on alpha" << std::endl;
            double m = settingsPtr->getPriorConcMean();
            double v = settingsPtr->getPriorConcVariance();
            gammaBeta  = m / v;
            gammaAlpha = gammaBeta * m;
            alpha = Probability::Gamma::rv(&rng, gammaAlpha, gammaBeta);
            std::cout << "   *    Initial Alpha = " << alpha << std::endl;
            }
        
        for (int i=0; i<numPatrons; i++)
            {
            double newTableProb = alpha / (alpha + i);
            double u1 = rng.uniformRv();
            if (u1 < newTableProb)
                {
                // add a new table
                Table* t = addTable();
                t->addPatron(i);
                }
            else
                {
                // add the patron to one of the existing tables
                double u2 = rng.uniformRv();
                double sum = 0.0;
                for (Table* t : tables)
                    {
                    sum += (double)t->getNumPatrons() / i;
                    if (u2 < sum)
                        {
                        t->addPatron(i);
                        break;
                        }
                    }
                }
            }
        std::cout << "   *    Initial RGF = " << rgf() << std::endl;

        }
}

Restaurant::~Restaurant(void) {

    delete parameter;
}

Table* Restaurant::addTable(void) {

    Table* t = TableFactory::tableFactory().getTable();
    Parameter* newParm = parameter->newRandomInstance();
    newParm->setTable(t);
    t->addParameter( newParm );
    t->setRestaurant(this);
    tables.insert(t);
    return t;
}

Table* Restaurant::addAuxiliaryTable(void) {

    Table* t = TableFactory::tableFactory().getTable();
    Parameter* newParm = parameter->newRandomInstance();
    newParm->setTable(t);
    t->addParameter( newParm );
    t->setRestaurant(this);
    return t;
}

double Restaurant::calculateAlphaFromExpectedNumberOfTables(double expT, int np) {

    if (expT > np)
        Msg::error("The expected number of tables cannot be larger than the number of patrons (" + std::to_string(np) + "<" + std::to_string(expT) + ")");
     if (expT <= 1.0)
        Msg::error("The expected number of tables cannot be less than one");
       
    double a = 0.000001;
    double ea = expectedNumberOfTables(a, np);
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
        ea = expectedNumberOfTables(a, np);
        }
    //cout << ea << " <-> " << expT << " " << "alpha=" << a << endl;
    return a;
}

Table* Restaurant::chooseTable(std::map<Table*,double>& lnProbs) {

    RandomVariable& rv = RandomVariable::randomVariableInstance();
    double u = rv.uniformRv();
    double sum = 0.0;
    for (std::map<Table*,double>::iterator it = lnProbs.begin(); it != lnProbs.end(); it++)
        {
        sum += it->second;
        if (u < sum)
            return it->first;
        }
    return NULL;
}

Parameter* Restaurant::copyParameter(Parameter* parmToCopy) {

    Parameter* parm = NULL;
    if (dynamic_cast<ParameterTree*>(parmToCopy) != NULL)
        parm = new ParameterTree(*dynamic_cast<ParameterTree*>(parmToCopy));
    else if (dynamic_cast<ParameterTreeLength*>(parmToCopy) != NULL)
        parm = new ParameterTreeLength(*dynamic_cast<ParameterTreeLength*>(parmToCopy));
    else if (dynamic_cast<ParameterBaseFrequencies*>(parmToCopy) != NULL)
        parm = new ParameterBaseFrequencies(*dynamic_cast<ParameterBaseFrequencies*>(parmToCopy));
    else if (dynamic_cast<ParameterExchangabilityRates*>(parmToCopy) != NULL)
        parm = new ParameterExchangabilityRates(*dynamic_cast<ParameterExchangabilityRates*>(parmToCopy));
    else if (dynamic_cast<ParameterGammaShape*>(parmToCopy) != NULL)
        parm = new ParameterGammaShape(*dynamic_cast<ParameterGammaShape*>(parmToCopy));
    return parm;
}

double Restaurant::expectedNumberOfTables(double a, int np) {

    double expectedNum = 0.0;
    for (int i=1; i<=np; i++)
        expectedNum += ( 1.0 / (i - 1.0 + a) );
    expectedNum *= a;
    return expectedNum;
}

Table* Restaurant::findTableWithPatron(int idx) {

    for (Table* t : tables)
        {
        if (t->hasPatron(idx) == true)
            return t;
        }
    return NULL;
}

void Restaurant::normalize(std::map<Table*,double>& lnProbs) {

    // find maximum value
    bool first = true;
    double maxVal = 0.0;
    for (std::map<Table*,double>::iterator it = lnProbs.begin(); it != lnProbs.end(); it++)
        {
        if (first == true)
            {
            maxVal = it->second;
            first = false;
            }
        else
            {
            if (it->second > maxVal)
                maxVal = it->second;
            }
        }
        
    // rescale
    for (std::map<Table*,double>::iterator it = lnProbs.begin(); it != lnProbs.end(); it++)
        it->second -= maxVal;
        
    // exponentiate
    double sum = 0.0;
    for (std::map<Table*,double>::iterator it = lnProbs.begin(); it != lnProbs.end(); it++)
        {
        double lnP = it->second;
        double p = exp(lnP);
        it->second = p;
        sum += p;
        }
        
    // rescale to probs
    for (std::map<Table*,double>::iterator it = lnProbs.begin(); it != lnProbs.end(); it++)
        {
        double p = it->second / sum;
        it->second = p;
        }

#   if 0
    for (std::map<Table*,double>::iterator it = lnProbs.begin(); it != lnProbs.end(); it++)
        std::cout << it->first << " " << it->second << std::endl;
#   endif
}

void Restaurant::print(void) {

    std::cout << "Restaurant " << parameter->getName() << " (" << rgf() << ")" << std::endl;
    int i = 0;
    for (Table* t : tables)
        {
        Parameter* parm = t->getParameter();
        std::cout << "   Table " << i << ": ";
        std::cout << parm->getValuesAsString(5);
        i++;

        std::cout << " [ ";
        std::set<int>& patrons = t->getPatrons();
        for (int p : patrons)
            std::cout << p << " ";
        std::cout << "]" << std::endl;
        }
}

std::string Restaurant::rgf(void) {

    std::vector<int> rgf(numPatrons);
    for (int i=0; i<numPatrons; i++)
        rgf[i] = 0;
    int idx = 0;
    for (int i=0; i<numPatrons; i++)
        {
        if (rgf[i] == 0)
            {
            idx++;
            Table* t = findTableWithPatron(i);
            if (t == NULL)
                Msg::error("Could not find patron " + std::to_string(i) + " in the restaurant (" + parameter->getName() + ")");
            std::set<int> tPatrons = t->getPatrons();
            for (int p : tPatrons)
                {
                rgf[p] = idx;
                }
            }
        }
    
    std::string rgfStr = "";
    for (int i=0; i<numPatrons; i++)
        {
        if (i != 0)
            rgfStr += ",";
        rgfStr += std::to_string(rgf[i]);
        }
    return rgfStr;
}

void Restaurant::removeTable(Table* tab) {

    TableFactory& tf = TableFactory::tableFactory();
    tf.returnToPool(tab);
    tables.erase(tab);
}

double Restaurant::sampleAlpha(int k, int n, double oldAlpha, double a, double b) {

    RandomVariable& rng = RandomVariable::randomVariableInstance();

	/* Step 1: Draw a Beta(oldAlpha+1, n) distribution to get eta */
	std::vector<double> z(2);
	std::vector<double> f(2);
	z[0] = oldAlpha + 1.0;
	z[1] = (double)n;
    Probability::Dirichlet::rv(&rng, z, f);
	double eta = f[0];
	
	/* Step 2: Draw a new value for alpha, based on k and eta */
	double u = rng.uniformRv();
	double x = ( a + (double)k - 1.0 ) / ( (double)n * (b - log(eta)) );
	double newAlpha;
	if ( (u / (1.0 - u)) < x)
		newAlpha = Probability::Gamma::rv(&rng, a + k, b - log(eta));
	else
		newAlpha = Probability::Gamma::rv(&rng, a + k - 1.0, b - log(eta));
		
	return newAlpha;
}

double Restaurant::update(void) {

    // The update of a restaurant does two thigns. First, we visit all of the
    // current tables in the restaurant, and update the parameter on each. We
    // then reseat patrons at tables according to Algorithm 8 of Neal (2000).
        
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    
    // calculate the initial likelihoods
    std::vector<double> curLnL(numPatrons);
    std::vector<double> newLnL(numPatrons);
    for (int i=0; i<numPatrons; i++)
        {
        Chunk* chunk = modelPtr->getChunk(i);
        chunk->updateRateMatrix();
        chunk->updateTransitionProbabilities();
        curLnL[i] = chunk->lnLikelihood();
        }
        
    // 1. update the parameters on each table
    for (Table* tab : tables)
        {
        // change the parameter
        Parameter* parm = tab->getParameter();
        double lnPriorRatio = -parm->lnProbability();
        double lnProposalRatio = parm->update();
        lnPriorRatio += parm->lnProbability();
        
        // calculate likelihood
        std::set<int>& tablePatrons = tab->getPatrons();
        double lnLikelihoodRatio = 0.0;
        for (int i : tablePatrons)
            {
            Chunk* chunk = modelPtr->getChunk(i);
            if (parm->getUpdateModifiesEigens() == true)
                {
                chunk->flipActiveEigens();
                chunk->updateRateMatrix();
                }
            chunk->flipActiveTransitionProbabilities();
            chunk->updateTransitionProbabilities();

            newLnL[i] = modelPtr->lnLikelihood(i);
            lnLikelihoodRatio += (newLnL[i] - curLnL[i]);
            }
            
        // accept/reject
        double lnR = lnLikelihoodRatio + lnPriorRatio + lnProposalRatio;
        bool accept = false;
        if (log(rng.uniformRv()) < lnR)
            accept = true;
            
        // update
        if (accept == true)
            {
            // proposed update is accepted
            parm->accept();
            for (int i : tablePatrons)
                curLnL[i] = newLnL[i];
            UpdateInfo::updateInfo().accept();
            }
        else
            {
            // proposed update is rejected
            parm->reject();
            for (int i : tablePatrons)
                {
                Chunk* chunk = modelPtr->getChunk(i);
                if (parm->getUpdateModifiesEigens() == true)
                    chunk->flipActiveEigens();
                chunk->flipActiveTransitionProbabilities();
                }
            UpdateInfo::updateInfo().reject();
            }
        }

    // 2. loop over all of the patrons, reseating each (alg. 8)
    if (isSeatingRv == true)
        {
        std::set<Table*> auxiliaryTables;
        int numAuxliaryTables = 10; // was 3
        for (int n=0; n<numPatrons; n++)
            {
            // get the chunk for the patron
            Chunk* chunk = modelPtr->getChunk(n);
            
            // remove the patron from its current table
            Table* tab = findTableWithPatron(n);
            if (tab == NULL)
                Msg::error("Could not find table with patron " + std::to_string(n));
            tab->removePatron(n);
            
            // if there are no longer any patrons at the table, remove the table from the restaurant
            if (tab->getNumPatrons() == 0)
                removeTable(tab);
            
            // make numAuxiliary tables drawing the parameter values from the prior distribution
            // the auxiliary tables are added to the tables set
            for (int i=0; i<numAuxliaryTables; i++)
                auxiliaryTables.insert( addTable() );
                            
            // calculate the likelihood when the patron is seated at each of the tables (including the auxiliary tables)
            std::map<Table*,double> tableLikes;
            for (Table* t : tables)
                {
                t->addPatron(n);
                if (t->getParameter()->getUpdateModifiesEigens() == true)
                    chunk->updateRateMatrix();
                chunk->updateTransitionProbabilities();
                double lnL = chunk->lnLikelihood();
                std::set<Table*>::iterator it = auxiliaryTables.find(t);
                double lnP = 0.0;
                if (it == auxiliaryTables.end())
                    lnP = log( (double)(t->getNumPatrons() - 1) ); // subtract the patron that was added for calculation purposes
                else
                    lnP = log(alpha/numAuxliaryTables);
                tableLikes.insert( std::make_pair(t, lnL+lnP) );
                t->removePatron(n);
                }

            // choose a table in proportion to the likelihood to reseat
            normalize(tableLikes);
            tab = chooseTable(tableLikes);
            if (tab == NULL)
                Msg::error("Could not find a new table for reseating patron " + std::to_string(n));
            tab->addPatron(n);
            
            // sort out the likelihood for the patron
            if (tab->getParameter()->getUpdateModifiesEigens() == true)
                chunk->updateRateMatrix();
            chunk->updateTransitionProbabilities();
            
            // remove any unused tables
            std::vector<Table*> emptyTables;
            for (Table* t : tables)
                {
                if (t->getNumPatrons() == 0)
                    emptyTables.push_back(t);
                }
            for (int i=0; i<emptyTables.size(); i++)
                removeTable(emptyTables[i]);
            auxiliaryTables.clear();
            }
        }
 
    // 3. refresh all likelihood calculations (necessary?)
    double lnL = 0.0;
    for (int i=0; i<numPatrons; i++)
        {
        Chunk* chunk = modelPtr->getChunk(i);
        chunk->updateRateMatrix();
        chunk->updateTransitionProbabilities();
        lnL += chunk->lnLikelihood();
        }
        
    // 4. update the concentration parameter (Gibb's sampler)
	if ( settingsPtr->getIsConcentrationParameterFixed() == false && isSeatingRv == true )
		alpha = sampleAlpha( getNumTables(), numPatrons, alpha, gammaAlpha, gammaBeta );

    return lnL;
}
