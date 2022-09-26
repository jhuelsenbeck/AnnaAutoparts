#ifndef SubRates_H
#define SubRates_H

#include "Parm.h"
#include <string>



class MbRandom;
class SubRates : public Parm {

	public:
                	              SubRates(MbRandom *rp, Model *mp, std::string nm, double tn);
								  SubRates(SubRates &s);
                	              ~SubRates(void);
					   SubRates   &operator=(SubRates &s);
					     double   update(void);
						 double   lnPriorProb(void);
					std::string   getParmString(int n);
					       void   clone(SubRates &s);
				           void   print(void);
						 double   getSubRate(int i) { return r[i]; }
		   std::vector<double>&   getSubRate(void) { return r; }
					std::string   getParmHeader(int n);
                	              
	private:
						 double   normalizeRates(std::vector<double> &a, double minVal, double total);
	                     double   alpha0;
			std::vector<double>   a;
			std::vector<double>   r;
                            
};

#endif