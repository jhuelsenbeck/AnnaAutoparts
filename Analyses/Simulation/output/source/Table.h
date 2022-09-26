#ifndef Table_H
#define Table_H

#include <set>
#include <string>

enum ParmId { ASRV, TREE, SUBRATE, BASEFREQ, LENGTH };

class Alignment;
class MbRandom;
class Model;
class Parm;
class Settings;

class Table {

	public:
                            Table(MbRandom* rp, Settings* sp, Alignment* ap, Model* mp, ParmId pid);
						   ~Table(void);
					 void   flipCurParmId(void) { (curParmId == 0 ? curParmId = 1 : curParmId = 0); }
					Parm*   getParm(void) { return parameter[curParmId]; }
			  std::string   getParmString(int n);
		   std::set<int>&   getPatrons(void) { return patrons; }
					 bool   isPatronAtTable(int i);
					  int   numPatronsAtTable(void) { return patrons.size(); }
					 void   print(void);
					 void   removePatron(int i);
					 void   seatPatron(int i);

	private:
	                 void   makeNewParm(void);
	        std::set<int>   patrons;
			   Alignment*   alignmentPtr;
				MbRandom*   ranPtr;
				   Model*   modelPtr;
				Settings*   settingsPtr;
					Parm*   parameter[2];
				      int   curParmId;
				   ParmId   parmId;
};


#endif