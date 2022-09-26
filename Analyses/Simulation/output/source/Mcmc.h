#ifndef Mcmc_H
#define Mcmc_H

#include <valarray>
#include <fstream>

class MbRandom;
class Model;
class Settings;

class Mcmc {

	public:
                	              Mcmc(Settings* sp, Model* mp, MbRandom* rp);  
								 ~Mcmc(void);
                	              
	private:
				         double   aicm(void);
				         double   marginalLikelihood(void);
                         double   Pr_harmonic(const std::valarray<double>& v);
                         double   Pr_smoothed(const std::valarray<double>& v);
                         double   Pr_smoothed(const std::valarray<double>& v, double delta, double Pdata);
	                       void   saveStates(int n, double lnL);
						 Model*   modelPtr;
					  Settings*   settingsPtr;
					  MbRandom*   ranPtr;
				  std::ofstream   lnOut;
			std::vector<double>   logLikes;
                            
};

#endif