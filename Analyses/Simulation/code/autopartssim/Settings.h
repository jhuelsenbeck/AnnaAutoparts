#ifndef Settings_H
#define Settings_H

#include <string>
#include <vector>



class Settings {

	public:
                                        Settings(std::string pf);
								  int   getNumReplicates(void) { return numReplicates; }
                                  int   getNumTaxa(void) { return numTaxa; }
                    std::vector<int>&   getNumSites(void) { return numSites; }
                 std::vector<double>&   getBranchLengthPrior(void) { return branchLengthPrior; }
                 std::vector<double>&   getAlpha(void) { return alpha; }
  std::vector< std::vector<double> >&   getBaseFreqs(void) { return baseFreqs; }
  std::vector< std::vector<double> >&   getSubstitutionRates(void) { return substitutionRates; }

	private:
	                              int   numReplicates;
                                  int   numTaxa;
                     std::vector<int>   numSites;
                  std::vector<double>   branchLengthPrior;
                  std::vector<double>   alpha;
   std::vector< std::vector<double> >   baseFreqs;
   std::vector< std::vector<double> >   substitutionRates;
};

#endif