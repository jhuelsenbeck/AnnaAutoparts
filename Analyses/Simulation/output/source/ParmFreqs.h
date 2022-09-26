#ifndef BaseFreqs_H
#define BaseFreqs_H

#include "Parm.h"
#include <string>
#include <vector>



class MbRandom;
class BaseFreqs : public Parm {

	public:
                	              BaseFreqs(MbRandom *rp, Model *mp, std::string nm, double tn);
								  BaseFreqs(BaseFreqs &b);
								 ~BaseFreqs(void);
					  BaseFreqs   &operator=(BaseFreqs &b);
					     double   update(void);
					       void   clone(BaseFreqs &b);
						 double   lnPriorProb(void);
					std::string   getParmString(int n);
				           void   print(void);
						 double   getFreq(int i) { return f[i]; }
		   std::vector<double>&   getFreq(void) { return f; }
					std::string   getParmHeader(int n);
                	              
	private:
	                     double   alpha0;
			std::vector<double>   a;
			std::vector<double>   f;
                            
};

#endif