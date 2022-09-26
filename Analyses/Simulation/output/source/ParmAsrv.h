#ifndef ASRV_H
#define ASRV_H

#include "Parm.h"
#include <string>
#include <vector>



class MbRandom;
class Asrv : public Parm {

	public:
                	              Asrv(MbRandom *rp, Model *mp, std::string nm, double lm, int nc, double tn);  
								  Asrv(Asrv &a);
                	              ~Asrv(void);
						  Asrv&   operator=(Asrv &a);
						   void   clone(Asrv &a);
						 double   lnPriorProb(void);
					std::string   getParmString(int n);
						 double   update(void);
				           void   print(void);
						 double   getAlpha(void) { return alpha; }
						 double   getRate(int i) { return r[i]; }
		   std::vector<double>&   getRate(void) { return r; }
			                int   getNumGammaCats(void) { return numCats; }
					std::string   getParmHeader(int n);
                	              
	private:
					     double   alpha;
						 double   tuning;
						    int   numCats;
		    std::vector<double>   r;
						 double   lambda; // the exponential parameter for the prior on alpha
                            
};

#endif