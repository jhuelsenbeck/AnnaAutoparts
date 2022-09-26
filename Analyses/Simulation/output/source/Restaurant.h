#ifndef Restaurant_H
#define Restaurant_H

#include <fstream>
#include <map>
#include <set>
#include <vector>
#include "Table.h"

class Alignment;
class MbRandom;
class Model;
class Parm;
class Partition;
class Settings;

class Restaurant {

	public:
                            Restaurant(MbRandom* rp, Settings* sp, Alignment* ap, Model* mp, ParmId pid, int np, double ps);
						   ~Restaurant(void);
				     void   addPartition(void);
					 bool   change(void);
					 bool   changeParmOnTable(Table* tbl);
			  std::string   getName(void) { return name; }
			  std::string   getShortName(void) { return shortName; }
					  int   getNumTables(void) { return tables.size(); }
				   ParmId   getParmId(void) { return parmId; }
				   Table*   getTableWithPatron(int patron);
					 void   print(void);
					 void   saveState(int n);
					 void   simulateFromPrior(std::vector<int>& x, int nReps);
					 void   summarizePartitions(void);
					 void   updateSeating(int patron);

	private:
	               double   acceptanceProb(double lnR);
				   double   calcAlphaFromEt(double expT);
				     void   deleteUnoccupiedTables(void);
				   double   expNumTables(double a);
				     void   normalizeLogProbabilitiesInMap(std::map<Table*,double>& m);
                     void   mySort(std::vector<int> &item, std::vector<int> &assoc, int count);
                     void   sort2(std::vector<int> &item, std::vector<int> &assoc, int left, int right);
				   Table*   pickTableAtRandomFromPrior(void);
				   Table*   pickTableUniformlyAtRandom(void);
				     void   removeTable(Table* tbl);
				   double   sampleAlpha(int k, int n, double oldAlpha, double a, double b);
					 void   setRgf(void);
		 std::set<Table*>   tables;
			   Alignment*   alignmentPtr;
				MbRandom*   ranPtr;
				   Model*   modelPtr;
				Settings*   settingsPtr;
				   ParmId   parmId;
			  std::string   name;
			  std::string   shortName;
                      int   numPatrons;
				   double   alpha;
				   double   gammaAlpha;
				   double   gammaBeta;
				   double   probUpdateSeating;
					  int   numAuxTables;
			std::ofstream   parmOut;
			std::ofstream   partOut;
					 int*   rgf;
  std::vector<Partition*>   sampledPartitions;
};


#endif