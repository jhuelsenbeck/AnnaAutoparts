#include "MbRandom.h"
#include "ParmSubrates.h"
#include <iostream>

#define MIN_RATE	0.0000001

#undef DEBUG_SUBRATES




SubRates::SubRates(MbRandom *rp, Model *mp, std::string nm, double tn) : Parm(rp, mp, nm) {

	alpha0 = tn;
	r = std::vector<double>( 6 );
	a = std::vector<double>( 6 );
	for (int i=0; i<6; i++)
		a[i] = 1.0;
	ranPtr->dirichletRv(a, r);
	normalizeRates(r, MIN_RATE, 1.0);	
}

SubRates::SubRates(SubRates &s) : Parm(s.ranPtr, s.modelPtr, s.parmName) {

	r = std::vector<double>( 6 );
	a = std::vector<double>( 6 );
	clone( s );

}

SubRates::~SubRates(void) {

}

SubRates &SubRates::operator=(SubRates &s) {

	if (this != &s)
		{
		Parm::operator=(s);
		}
	return *this;

}

double SubRates::lnPriorProb(void) {

	return ranPtr->lnDirichletPdf(a, r);
	
}

std::string SubRates::getParmString(int n) {

	char temp[200];
	sprintf(temp, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t", r[0], r[1], r[2], r[3], r[4], r[5]);
	std::string tempString = temp;
	return tempString;
	
}

std::string SubRates::getParmHeader(int n) {

	char temp[500];
	if  ( n == -1 )
		sprintf (temp, "R[AC]\tR[AG]\tR[AT]\tR[CG]\tR[CT]\tR[GT]\t");
	else
		sprintf (temp, "R[AC,%d]\tR[AG,%d]\tR[AT,%d]\tR[CG,%d]\tR[CT,%d]\tR[GT,%d]\t", n, n, n, n, n, n);
	std::string tempString = temp;
	return tempString;
}

double SubRates::update(void) {

#if 0

#	if defined (DEBUG_SUBRATES)
	cout << r << " -> ";
#	endif
	std::vector<double> aForward(6);
	std::vector<double> aReverse(6);
	std::vector<double> oldRates(6);
	for (int i=0; i<6; i++)
		{
		oldRates[i] = r[i];
		aForward[i] = r[i] * alpha0;
		}
	ranPtr->dirichletRv(aForward, r);
	for (int i=0; i<6; i++)
		aReverse[i] = r[i] * alpha0;
#	if defined (DEBUG_SUBRATES)
	cout << r << endl;
#	endif
	return ranPtr->lnDirichletPdf(aReverse, oldRates) - ranPtr->lnDirichletPdf(aForward, r);
	
#else

#	if defined (DEBUG_SUBRATES)
	cout << r << " -> ";;
#	endif

	/* how many substitution rates do we update? */
	int k = 2;

	/* pick the substitution rates to update */
	std::vector<int> valueToUpdate(k);
	for (int i=0; i<k; i++)
		{
		
		int x;
		bool isThisOk;
		do
			{
			x = (int)( ranPtr->uniformRv() * 6 );
			isThisOk = true;
			for (int j=0; j<i; j++)
				{
				if ( x == valueToUpdate[j] )
					{
					isThisOk = false;
					break;
					}
				}
			} while( isThisOk == false );
		valueToUpdate[i] = x;
		}

	/* pick new proportions from a Dirichlet */
	std::vector<double> srNew(k+1);
	std::vector<double> srOld(k+1);
	std::vector<double> aNew(k+1);
	std::vector<double> aOld(k+1);
	double sum = 0.0;
	for (int i=0; i<k; i++) 
		{
		aNew[i] = r[ valueToUpdate[i] ] * alpha0;
		srOld[i] = r[ valueToUpdate[i] ];
		sum += r[ valueToUpdate[i] ];
		}
	aNew[k] = alpha0  * (1.0 - sum);
	srOld[k] = 1.0 - sum;
	ranPtr->dirichletRv(aNew, srNew);
	sum = normalizeRates(srNew, MIN_RATE, 1.0);
	
	for (int i=0; i<k+1; i++)
		aOld[i] = srNew[i] * alpha0;
		
	for (int id=0; id<6; id++)
		{
		bool isThisAnUpdatedRate = false;
		int updateElement = 0;
		for (int j=0; j<k; j++)
			{
			if ( id == valueToUpdate[j] )
				{
				isThisAnUpdatedRate = true;
				updateElement = j;
				break;
				}
			}
		if (isThisAnUpdatedRate == true)
			r[id] = srNew[ updateElement ];
		else
			r[id] *= srNew[k] / srOld[k];
		}

#	if defined (DEBUG_SUBRATES)
	cout << r << endl;
#	endif

	/* calculate prior ratio */
	double lnProposalProb = ( ranPtr->lnDirichletPdf(aOld, srOld) - ranPtr->lnDirichletPdf(aNew, srNew) ) + (6 - k + 1) * (log(srNew[k]) - log(srOld[k]) );
	return lnProposalProb;
	
#endif
}

void SubRates::clone(SubRates &s) {

	alpha0 = s.alpha0;
	for (int i=0; i<6; i++)
		{
		a[i] = s.a[i];
		r[i] = s.r[i];
		}

}

double SubRates::normalizeRates(std::vector<double> &a, double minVal, double total) {

	int n = a.size();
	double normalizeTo = total;
	double sum = 0.0;
	for (int i=0; i<n; i++)
		{
		if (a[i] <= minVal)
			normalizeTo -= minVal;
		else
			sum += a[i];
		}
	for (int i=0; i<n; i++)
		{
		if (a[i] <= minVal)
			a[i] = minVal;
		else
			a[i] *= (normalizeTo/sum);
		}
	sum = 0.0;
	for (int i=0; i<n; i++)
		sum += a[i];
	return sum;
}

void SubRates::print(void) {

	std::cout << "SubRates = ";
	for (int i=0; i<r.size(); i++)
		std::cout << r[i] << " ";
	std::cout << std::endl;
}
