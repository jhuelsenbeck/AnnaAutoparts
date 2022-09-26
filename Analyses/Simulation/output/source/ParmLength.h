#ifndef TreeLength_H
#define TreeLength_H

#include "Parm.h"

class MbRandom;
class TreeLength : public Parm {

	public:
                	              TreeLength(MbRandom *rp, Model *mp, std::string nm, double lm, int nb, double tn);  
								  TreeLength(TreeLength &t);
								 ~TreeLength(void);
					 TreeLength   &operator=(TreeLength &t);
						   void   clone(TreeLength &t);
						 double   lnPriorProb(void);
					std::string   getParmString(int n);
						 double   update(void);
				           void   print(void);
						 double   getLength(void) { return length; }
					std::string   getParmHeader(int n);
                	              
	private:
					     double   length;
						 double   tuning;
						    int   numBranches;
						 double   lambda; // the exponential parameter for the prior on branch lengths
                            
};

#endif