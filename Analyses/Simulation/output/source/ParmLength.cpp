#include <cmath>
#include <iostream>
#include "MbRandom.h"
#include "ParmLength.h"



TreeLength::TreeLength(MbRandom *rp, Model *mp, std::string nm, double lm, int nb, double tn) : Parm(rp, mp, nm) {

	numBranches = nb;
	lambda      = lm;
	tuning      = tn;
	length      = ranPtr->gammaRv( numBranches, lambda );
}

TreeLength::TreeLength(TreeLength &t) : Parm(t.ranPtr, t.modelPtr, t.parmName) {

	numBranches = t.numBranches;
	lambda      = t.lambda;
	tuning      = t.tuning;
	length      = t.length;
	
}

TreeLength::~TreeLength(void) {

}

TreeLength &TreeLength::operator=(TreeLength &t) {

	if (this != &t)
		{
		Parm::operator=(t);
		}
	return *this;

}

void TreeLength::clone(TreeLength &t) {

	numBranches = t.numBranches;
	lambda      = t.lambda;
	tuning      = t.tuning;
	length      = t.length;

}

std::string TreeLength::getParmString(int n) {

	char temp[100];
	sprintf(temp, "%lf\t", length);
	std::string tempString = temp;
	return tempString;
	
}

std::string TreeLength::getParmHeader(int n) {

	char temp[1000];
	if  ( n == -1 )
		sprintf (temp, "TL\t");
	else
		sprintf (temp, "TL[%d]\t", n);
	std::string tempString = temp;
	return tempString;
}

double TreeLength::lnPriorProb(void) {

	return ranPtr->lnGammaPdf( numBranches, lambda, length );
	
}

void TreeLength::print(void) {

	std::cout << "Tree Length = " << length << std::endl;
	
}

double TreeLength::update(void) {

	double oldValue = length;
	double newValue = oldValue * exp( tuning*(ranPtr->uniformRv()-0.5) );
	length = newValue;
	return log(newValue) - log(oldValue);
	
}


