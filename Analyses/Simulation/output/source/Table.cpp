#include "Alignment.h"
#include "ParmAsrv.h"
#include "ParmFreqs.h"
#include "ParmLength.h"
#include "ParmSubrates.h"
#include "ParmTree.h"
#include "Settings.h"
#include "Table.h"
#include "Util.h"



Table::Table(MbRandom* rp, Settings* sp, Alignment* ap, Model* mp, ParmId pid) {

	// remember the location of some objects
	ranPtr       = rp;
	settingsPtr  = sp;
	alignmentPtr = ap;
	modelPtr     = mp;
	parmId       = pid;
	curParmId    = 0;
	makeNewParm();}

Table::~Table(void) {

	if (parmId == ASRV)
		{
		Asrv *derivedPtr0 = dynamic_cast<Asrv *>(parameter[0]);
		Asrv *derivedPtr1 = dynamic_cast<Asrv *>(parameter[1]);
		if ( derivedPtr0 != 0 && derivedPtr1 != 0 )
			{
			delete derivedPtr0;
			delete derivedPtr1;
			}
		else
			{
			std::cout << "ERROR: Problem deleting Asrv parameter" << std::endl;
			exit(1);
			}
		}
	else if (parmId == TREE)
		{
		Tree *derivedPtr0 = dynamic_cast<Tree *>(parameter[0]);
		Tree *derivedPtr1 = dynamic_cast<Tree *>(parameter[1]);
		if ( derivedPtr0 != 0 && derivedPtr1 != 0 )
			{
			delete derivedPtr0;
			delete derivedPtr1;
			}
		else
			{
			std::cout << "ERROR: Problem deleting Tree parameter" << std::endl;
			exit(1);
			}
		}
	else if (parmId == SUBRATE)
		{
		SubRates *derivedPtr0 = dynamic_cast<SubRates *>(parameter[0]);
		SubRates *derivedPtr1 = dynamic_cast<SubRates *>(parameter[1]);
		if ( derivedPtr0 != 0 && derivedPtr1 != 0 )
			{
			delete derivedPtr0;
			delete derivedPtr1;
			}
		else
			{
			std::cout << "ERROR: Problem deleting SubRates parameter" << std::endl;
			exit(1);
			}
		}
	else if (parmId == BASEFREQ)
		{
		BaseFreqs *derivedPtr0 = dynamic_cast<BaseFreqs *>(parameter[0]);
		BaseFreqs *derivedPtr1 = dynamic_cast<BaseFreqs *>(parameter[1]);
		if ( derivedPtr0 != 0 && derivedPtr1 != 0 )
			{
			delete derivedPtr0;
			delete derivedPtr1;
			}
		else
			{
			std::cout << "ERROR: Problem deleting BaseFreqs parameter" << std::endl;
			exit(1);
			}
		}
	else if (parmId == LENGTH)
		{
		TreeLength *derivedPtr0 = dynamic_cast<TreeLength *>(parameter[0]);
		TreeLength *derivedPtr1 = dynamic_cast<TreeLength *>(parameter[1]);
		if ( derivedPtr0 != 0 && derivedPtr1 != 0 )
			{
			delete derivedPtr0;
			delete derivedPtr1;
			}
		else
			{
			std::cout << "ERROR: Problem deleting TreeLength parameter" << std::endl;
			exit(1);
			}
		}
	patrons.clear();
}

std::string Table::getParmString(int n) {

	return parameter[curParmId]->getParmString(n);
}

bool Table::isPatronAtTable(int i) {
	
	std::set<int>::iterator it = patrons.find( i );
	if (it != patrons.end())
		return true;
	return false;
}

void Table::makeNewParm(void) {

	if (parmId == ASRV)
		{
		parameter[0] = new Asrv( ranPtr, modelPtr, "Asrv", settingsPtr->getAsrvLambda(), settingsPtr->getNumGammaCats(), settingsPtr->getTuningParm("asrv") );
		parameter[1] = new Asrv( ranPtr, modelPtr, "Asrv", settingsPtr->getAsrvLambda(), settingsPtr->getNumGammaCats(), settingsPtr->getTuningParm("asrv") );
		(*parameter[1]) = (*parameter[0]);
		}
	else if (parmId == TREE)
		{
		/* get the constraints tree */
		std::string treeStr = "";
		if ( settingsPtr->getTreeFileName() != "" )
			{
			treeStr = getLineFromFile( settingsPtr->getTreeFileName(), 0 );
			parameter[0] = new Tree( ranPtr, modelPtr, "Tree", alignmentPtr, settingsPtr->getBrlenLambda(), settingsPtr->getTuningParm("tree"), treeStr );
			parameter[1] = new Tree( ranPtr, modelPtr, "Tree", alignmentPtr, settingsPtr->getBrlenLambda(), settingsPtr->getTuningParm("tree"), treeStr );
			}
		else
			{
			parameter[0] = new Tree( ranPtr, modelPtr, "Tree", alignmentPtr, settingsPtr->getBrlenLambda(), settingsPtr->getTuningParm("tree") );
			parameter[1] = new Tree( ranPtr, modelPtr, "Tree", alignmentPtr, settingsPtr->getBrlenLambda(), settingsPtr->getTuningParm("tree") );
			}
		(*parameter[1]) = (*parameter[0]);
		}
	else if (parmId == SUBRATE)
		{
		parameter[0] = new SubRates( ranPtr, modelPtr, "Rates", settingsPtr->getTuningParm("rates") );
		parameter[1] = new SubRates( ranPtr, modelPtr, "Rates", settingsPtr->getTuningParm("rates") );
		(*parameter[1]) = (*parameter[0]);
		}
	else if (parmId == BASEFREQ)
		{
		parameter[0] = new BaseFreqs( ranPtr, modelPtr, "Freqs", settingsPtr->getTuningParm("freqs") );
		parameter[1] = new BaseFreqs( ranPtr, modelPtr, "Freqs", settingsPtr->getTuningParm("freqs") );
		(*parameter[1]) = (*parameter[0]);
		}
	else if (parmId == LENGTH)
		{
		parameter[0] = new TreeLength( ranPtr, modelPtr, "Length", settingsPtr->getBrlenLambda(), 2*alignmentPtr->getNumTaxa()-3, settingsPtr->getTuningParm("length") );
		parameter[1] = new TreeLength( ranPtr, modelPtr, "Length", settingsPtr->getBrlenLambda(), 2*alignmentPtr->getNumTaxa()-3, settingsPtr->getTuningParm("length") );
		(*parameter[1]) = (*parameter[0]);
		}
}

void Table::print(void) {

	std::cout << "Patrons at table: \"";
	for (std::set<int>::iterator p = patrons.begin(); p != patrons.end(); p++)
		std::cout << (*p);
	std::cout << "\"" << std::endl;
	parameter[curParmId]->print();
}

void Table::removePatron(int i) {

	patrons.erase( i );
}

void Table::seatPatron(int i) {

	patrons.insert( i );
}

