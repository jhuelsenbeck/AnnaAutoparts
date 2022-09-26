#ifndef Model_H
#define Model_H

#include <vector>

class Alignment;
class Asrv;
class BaseFreqs;
class Chunk;
class MbRandom;
class Restaurant;
class Settings;
class StateSets;
class SubRates;
class Table;
class Tree;
class TreeLength;

class Model {

	public:
                            Model(Settings* sp, MbRandom* rp, Alignment* ap);
						   ~Model(void);
					Asrv*   findAsrv(int part);
			   BaseFreqs*   findBaseFreqs(int part);
			    SubRates*   findSubRates(int part);
					Tree*   findTree(int part);
			  TreeLength*   findTreeLength(int part);
                   Chunk*   getChunk(int i) { return chunks[i]; }
			          int   getNumRestaurants(void) { return restaurants.size(); }
			  Restaurant*   getRestaurant(int i) { return restaurants[i]; }
					  int   getRestaurantId(Restaurant* rest);
			  Restaurant*   getRestaurantToChange(void);
			   StateSets*   getStateSetPtr(void) { return stateSets; }
                   double   lnLikelihood(void);
			       double   lnLikelihood(int patron, bool storeScore);
				   double   lnLikelihood(bool storeScore);
                     void   setStoredLnLToMostRecentLnL(Table* tbl);
                     void   setUpdateFlags(Table* tbl);

	private:
				MbRandom*   ranPtr;
				Settings*   settingsPtr;
			   Alignment*   alignmentPtr;
	  std::vector<Chunk*>   chunks;
			   StateSets*   stateSets;
 std::vector<Restaurant*>   restaurants;
	  std::vector<double>   proposalProb;
	                  int   numSubsets;
};


#endif