#include "MbRandom.h"
#include "ParmFreqs.h"
#include <iostream>



#pragma mark Constructors

BaseFreqs::BaseFreqs(MbRandom *rp, Model *mp, std::string nm, double tn) : Parm(rp, mp, nm) {

	f = std::vector<double>( 4 );
	a = std::vector<double>( 4 );
	alpha0 = tn;
	for (int i=0; i<4; i++)
		a[i] = 1.0;
	ranPtr->dirichletRv(a, f);
		
}

BaseFreqs::BaseFreqs(BaseFreqs &b) : Parm(b.ranPtr, b.modelPtr, b.parmName) {

	f = std::vector<double>( 4 );
	a = std::vector<double>( 4 );
	clone( b );
		
}

#pragma mark Destructor

BaseFreqs::~BaseFreqs(void) {

}

#pragma mark Operators

BaseFreqs &BaseFreqs::operator=(BaseFreqs &b) {

	if (this != &b)
		{
		Parm::operator=(b);
		}
	return *this;

}

#pragma mark Class Functions

double BaseFreqs::lnPriorProb(void) {

	return ranPtr->lnDirichletPdf(a, f);
	
}

std::string BaseFreqs::getParmString(int n) {

	char temp[200];
	sprintf(temp, "%lf\t%lf\t%lf\t%lf\t", f[0], f[1], f[2], f[3]);
	std::string tempString = temp;
	return tempString;
	
}

std::string BaseFreqs::getParmHeader(int n) {

	char temp[500];
	if  ( n == -1 )
		sprintf (temp, "Pi[A]\tPi[C]\tPi[G]\tPi[T]\t");
	else
		sprintf (temp, "Pi[A,%d]\tPi[C,%d]\tPi[G,%d]\tPi[T,%d]\t", n, n, n, n);
	std::string tempString = temp;
	return tempString;
}

double BaseFreqs::update(void) {

	std::vector<double> aForward(4);
	std::vector<double> aReverse(4);
	std::vector<double> oldFreqs(4);
	for (int i=0; i<4; i++)
		{
		oldFreqs[i] = f[i];
		aForward[i] = f[i] * alpha0;
		}
	ranPtr->dirichletRv(aForward, f);
	for (int i=0; i<4; i++)
		aReverse[i] = f[i] * alpha0;
	return ranPtr->lnDirichletPdf(aReverse, oldFreqs) - ranPtr->lnDirichletPdf(aForward, f);
	
}

void BaseFreqs::clone(BaseFreqs &b) {

	alpha0 = b.alpha0;
	for (int i=0; i<4; i++)
		{
		a[i] = b.a[i];
		f[i] = b.f[i];
		}

}

void BaseFreqs::print(void) {

	std::cout << "Base Freqs = ";
	for (int i=0; i<f.size(); i++)
		std::cout << f[i] << " ";
	std::cout << std::endl;	
}

