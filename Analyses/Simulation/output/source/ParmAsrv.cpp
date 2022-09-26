#include <cmath>
#include <iostream>
#include "MbRandom.h"
#include "Parm.h"
#include "ParmAsrv.h"




Asrv::Asrv(MbRandom *rp, Model *mp, std::string nm, double lm, int nc, double tn) : Parm(rp, mp, nm) {

	lambda  = lm;
	tuning  = tn;
	alpha   = ranPtr->exponentialRv( lambda );
	numCats = nc;
	r       = std::vector<double>( numCats );
	ranPtr->discretizeGamma(r, alpha, alpha, numCats, false);
	
}

Asrv::Asrv(Asrv &a) : Parm(a.ranPtr, a.modelPtr, a.parmName) {

	r = std::vector<double>( a.numCats );
	clone(a);
		
}

Asrv::~Asrv(void) {

}

Asrv &Asrv::operator=(Asrv &a) {

	if (this != &a)
		{
		Parm::operator=(a);
		}
	return *this;

}

void Asrv::clone(Asrv &a) {

	lambda = a.lambda;
	tuning = a.tuning;
	alpha  = a.alpha;
	numCats = a.numCats;
	for (int i=0; i<numCats; i++)
		r[i] = a.r[i];
		
}

std::string Asrv::getParmString(int n) {

	char temp[100];
	sprintf(temp, "%lf\t", alpha);
	std::string tempString = temp;
	return tempString;
	
}

std::string Asrv::getParmHeader(int n) {

	char temp[1000];
	if  ( n == -1 )
		sprintf (temp, "Alpha\t");
	else
		sprintf (temp, "Alpha[%d]\t", n);
	std::string tempString = temp;
	return tempString;
}

void Asrv::print(void) {

	std::cout << "alpha = " << alpha << std::endl;
	
}

double Asrv::lnPriorProb(void) {

	return ranPtr->lnExponentialPdf( lambda, alpha );
	
}

double Asrv::update(void) {

	double oldValue = alpha;
	double newValue = oldValue * exp( tuning*(ranPtr->uniformRv()-0.5) );
	alpha = newValue;
	ranPtr->discretizeGamma(r, alpha, alpha, numCats, false);
	return log(newValue) - log(oldValue);
	
}
