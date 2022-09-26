#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <vector>
#include <valarray>
#include <algorithm>
#include <cmath>
#include "logsum.h"
#include "Mcmc.h"
#include "Model.h"
#include "Restaurant.h"
#include "Settings.h"



Mcmc::Mcmc(Settings* sp, Model* mp, MbRandom* rp) {

	// remember some pointers that we will need
	settingsPtr = sp;
	modelPtr    = mp;
	ranPtr      = rp;
	
	/* Open a file for printing the log likelihood. The parameters are saved through
	   the Restaurant objects. */
	std::string lnFileName = settingsPtr->getOutPutFileName() + ".lnL";
	lnOut.open( lnFileName.c_str(), std::ios::out );
	
    // allocate information for the acceptance rates
	int *numAccepted = new int[2 * modelPtr->getNumRestaurants()];
	int *numAttempted = numAccepted + modelPtr->getNumRestaurants();
	for (int i=0; i<modelPtr->getNumRestaurants(); i++)
		numAccepted[i] = numAttempted[i] = 0;
        
    // initialize the log likelihoods
    double oldLnL = modelPtr->lnLikelihood(true);
	
    // run the Markov chain
	std::cout << std::endl;
	std::cout << std::setw(6) << "Gen" << " -- lnL" << std::setw(11) << std::fixed << " --" << " Number of occupied tables" << std::setw(6) << "--" << std::fixed << " Parameter update" << std::fixed << std::endl;
	std::vector<double> mean(modelPtr->getNumRestaurants(), 0.0);
	for (int n=1; n<=settingsPtr->getChainLength(); n++)
		{
		// pick a parameter to change
		Restaurant* rest = modelPtr->getRestaurantToChange();
		int restaurantId = modelPtr->getRestaurantId(rest);
		
		// update the parameter
		bool wasAccepted = rest->change();
        
        double newLnL = modelPtr->lnLikelihood();
		
		if ( n >= settingsPtr->getBurnIn() )
			logLikes.push_back( newLnL );
		
		// print information to the screen
		if ( n % settingsPtr->getPrintFrequency() == 0 )
			{
			std::cout << std::setw(6) << n << " -- " << std::fixed << std::setprecision(3) << newLnL << " -- ";
			for (int i=0; i<modelPtr->getNumRestaurants(); i++)
				std::cout << modelPtr->getRestaurant(i)->getShortName() << "(" << modelPtr->getRestaurant(i)->getNumTables() << ") ";
			std::cout << "-- updating " << rest->getName();
			std::cout << std::endl;
			}
		for (int i=0; i<modelPtr->getNumRestaurants(); i++)
			mean[i] += modelPtr->getRestaurant(i)->getNumTables();
			
		// remember the state of the chain
		if ( n % settingsPtr->getSampleFrequency() == 0)
			saveStates( n, newLnL );
			
		// update information on acceptances
		numAttempted[restaurantId]++;
		if (wasAccepted == true)
			numAccepted[restaurantId]++;
			
		// sample the partitions for all of the restaurants
		if ( n >= settingsPtr->getBurnIn() )
			{
			for (int i=0; i<modelPtr->getNumRestaurants(); i++)
				modelPtr->getRestaurant(i)->addPartition();
			}
			
		oldLnL = newLnL;
		}
	
	for (int i=0; i<modelPtr->getNumRestaurants(); i++)
		std::cout << i << " -- " << std::fixed << std::setprecision(4) << mean[i] / settingsPtr->getChainLength() << std::endl;
	
	// print information on acceptances
	int longestName = 0;
	for (int i=0; i<modelPtr->getNumRestaurants(); i++)
		{
		std::string restaurantName = modelPtr->getRestaurant(i)->getName();
		if (restaurantName.size() > longestName)
			longestName = restaurantName.size();
		}
	std::cout << std::endl;
	std::cout << "Acceptance Rates:" << std::endl;
	std::cout << "Parameter";
	for (int j=0; j<longestName - 6; j++)
		std::cout << " ";
	std::cout << "  ";
	std::cout << std::setw(5) << "Tries" << "  ";
	std::cout << std::setw(5) << "Accep" << " ";
	std::cout << std::setw(5) << "Rate" << std::endl;
	for (int i=0; i<modelPtr->getNumRestaurants(); i++)
		{
		std::string restaurantName = modelPtr->getRestaurant(i)->getName();
		std::cout << restaurantName;
		for (int j=0; j<longestName - restaurantName.size(); j++)
			std::cout << " ";
		std::cout << "  ";
		std::cout << std::setw(5) << numAttempted[i] << "  ";
		std::cout << std::setw(5) << numAccepted[i] << "    ";
		if (numAttempted[i] > 0)
			std::cout << std::fixed << std::setprecision(2) << std::setw(5) << (double)numAccepted[i] / numAttempted[i];
		else
			std::cout << std::setw(5) << "N/A";
		std::cout << std::endl;
		}

	std::valarray<double> values( logLikes.size() );
	for (int i=0; i<values.size(); i++)
	  values[i] = logLikes[i];
	double ml1 = aicm();
	double ml2 = marginalLikelihood();
	double ml3 = Pr_harmonic(values);
	double ml4 = Pr_smoothed(values);
	std::cout << std::endl;
	std::cout << "Marginal lnL Estimates:"					<< std::endl;
	std::cout << "AICM                           = " << ml1 << std::endl;
	std::cout << "Harmonic mean (old way)        = " << ml2 << std::endl;
	std::cout << "Harmonic mean (Suchard)        = " << ml3 << std::endl;
	std::cout << "Marginal likelihood (smoothed) = " << ml4 << std::endl;
	std::cout << std::endl;
	
	// have every restaurant calculate a summary of its sampled partitions
	for (int i=0; i<modelPtr->getNumRestaurants(); i++)
		{
		modelPtr->getRestaurant(i)->summarizePartitions();
		}


	delete [] numAccepted;
}

Mcmc::~Mcmc(void) {

	lnOut.close();
}

void Mcmc::saveStates(int n, double lnL) {

	if (n == 0)
		lnOut << "Gen" << '\t' << "lnL" << std::endl;
	lnOut << n << '\t' << lnL << std::endl;
	
	for (int i=0; i<modelPtr->getNumRestaurants(); i++)
		{
		modelPtr->getRestaurant(i)->saveState(n);
		}
		
}

double Mcmc::aicm(void) {

	// find the number of samples
	int n = logLikes.size();
	
	if ( n <= 2 )
		return 0.0;
	
	// find maximum log likelihood
	double maxLnL = 0.0;
	for (int i=0; i<n; i++)
		{
		if (logLikes[i] > maxLnL || i == 0)
			maxLnL = logLikes[i];
		}
		
	// calculate the mean and variance of the log likelihood values
	double a = 0.0, aOld = 0.0, s = 0.0;
	for (int i=0; i<n; i++)
		{
		double x = logLikes[i];
		if (i == 0)
			{
			a = x;
			s = 0.0;
			}
		else
			{
			aOld = a;
			a = aOld + (x - aOld) / (i+1);
			s = s + (x - aOld) * (x - aOld);
			}
		}
	double m = a;
	double v = s / (logLikes.size()-1);
	double likeMax = m + v;
	if (likeMax < maxLnL)
		likeMax = maxLnL;
	double aicm = 2.0 * (m - v);
		
#	if 0
	for (int i=0; i<n; i++)
		cout << fixed << setprecision(6) << logLikes[i] << endl;
	cout << "m       = " << m << endl;
	cout << "v       = " << v << endl;
	cout << "likeMax = " << likeMax << endl;
	cout << "aicm    = " << aicm << endl;
	getchar();
#	endif

	return aicm;
}

double Mcmc::marginalLikelihood(void) {

	// find the number of samples
	int n = logLikes.size();
	
	if ( n <= 2 )
		return 0.0;
	
	// find minimum log likelihood
	double lnC = 0.0;
	for (int i=0; i<n; i++)
		{
		if (logLikes[i] < lnC)
			lnC = logLikes[i];
		}
		
	// calculate the harmonic mean, approximating the marginal likelihood
	double sum = 0.0;
	for (int i=0; i<n; i++)
		sum += exp( -logLikes[i] + lnC );
	double lnInverseH = -lnC + log((double)1.0/n) + log(sum);
	//lnMarginalLikelihood = -lnInverseH;
	
#	if 0
	cout << -lnInverseH << " ";
	for (int i=burn; i<n; i++)
		cout << logLikes[i] << ",";
	cout << endl;
#	endif
	return -lnInverseH;
}

double Mcmc::Pr_harmonic(const std::valarray<double>& v) {

	double sum = 0.0;
	for(int i=0; i<v.size(); i++)
		sum += v[i];
	double denominator = log_0;
	for (int i=0; i<v.size(); i++)
		denominator = logsum(denominator, sum-v[i]);
	return sum - denominator + log( double(v.size()) );
}

double Mcmc::Pr_smoothed(const std::valarray<double>& v,double delta,double Pdata) {

	double log_delta     = log(delta);
	double log_inv_delta = log(1.0-delta);
	double n             = v.size();
	double top           = log(n) + log_delta - log_inv_delta + Pdata;
	double bottom        = log(n) + log_delta - log_inv_delta;
	for(int i=0; i<v.size(); i++) 
		{
		double weight = -logsum(log_delta, log_inv_delta + v[i]-Pdata);
		top           = logsum(top, weight + v[i]);
		bottom        = logsum(bottom, weight);
		}
	return top - bottom;
}

double Mcmc::Pr_smoothed(const std::valarray<double>& v) {

	// sample from the prior w/ probability delta...
	double delta = 0.01;

	// use the harmonic estimator as the starting guess
	double Pdata = Pr_harmonic(v);

	// initialize this to a value which allows it to enter the loop
	double deltaP = 1.0;

	int iterations = 0;
	double dx = 10.0;
	while (std::abs(deltaP) > 1.0e-3) 
		{
		double g1 = Pr_smoothed(v, delta, Pdata) - Pdata;
		double Pdata2 = Pdata + g1;
		dx = g1 * 10;
		double g2 = Pr_smoothed(v, delta, Pdata+dx) - (Pdata+dx);
		double dgdx = (g2 - g1) / dx;

		double Pdata3 = Pdata - g1/dgdx;
		if ( Pdata3 < 2.0*Pdata || Pdata3 > 0 || Pdata3 > 0.5*Pdata )
			Pdata3 = Pdata + 10*g1;

		double g3 = Pr_smoothed(v, delta, Pdata3) - Pdata3;

		if ( std::abs(g3) <= std::abs(g2) && ((g3 > 0) || (std::abs(dgdx)>0.01)) ) 
			{
			// try to do Newton's method
			deltaP = Pdata3 - Pdata;
			Pdata = Pdata3;
			}
		else if ( std::abs(g2) <= std::abs(g1) ) 
			{
			// o/w see if we can go 10 times as far as one step
			Pdata2 += g2;
			deltaP = Pdata2 - Pdata;
			Pdata = Pdata2;
			}
		else 
			{
			// o/w go one step
			deltaP = g1;
			Pdata += g1;
			}

		iterations++;

		// if we aren't converging, warn and don't give an answer
		if (iterations > 400) 
			{
			std::cerr << "ERROR: Probabilities not converging!!!" << std::endl;
			return log_0;
			}

		}
	return Pdata;
}

