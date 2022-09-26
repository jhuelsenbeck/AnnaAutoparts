#include <cmath>
#include <iostream>
#include <iomanip>
#include "Alignment.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "StateSets.hpp"
#include "UserSettings.hpp"

#define UPPER  0
#define MIDDLE 1
#define LOWER  2



StateSets::StateSets(Alignment* ap, UserSettings* sp) {

	// set up pointers to other objects, and initialize some variables that will be commonly used
	alignmentPtr    = ap;
	settingsPtr     = sp;
	numTaxa         = alignmentPtr->getNumTaxa();
	numChar         = alignmentPtr->getNumSites();
	numNodes        = 2 * numTaxa - 2;
	oneStateSetSize = numNodes * numChar;
	
	// allocate the state sets vector
	sts = new unsigned[3 * oneStateSetSize];
	if ( !sts )
		{
		std::cerr << "ERROR: Problem allocating state sets" << std::endl;
		exit(1);
		}
	for (int i=0; i<3*oneStateSetSize; i++)
		sts[i] = 0;
		
	// and also allocate a matrix of pointers to elements of the sts vector, so we can
	// quickly position ourselves in that vector
	stsPtr = new unsigned**[3];
	stsPtr[0] = new unsigned*[3 * numNodes];
	stsPtr[1] = stsPtr[0] + numNodes;
	stsPtr[2] = stsPtr[1] + numNodes;
	for (int i=0; i<3; i++)
		{
		for (int j=0; j<numNodes; j++)
			{
			stsPtr[i][j] = &sts[ i*oneStateSetSize + j*(numChar) ];
			}
		}
		
	// and allocate a vector holding the lengths up to that point in the tree (state set)
	lengths = new int*[3];
	lengths[0] = new int[3 * numNodes];
	for (int i=1; i<3; i++)
		lengths[i] = lengths[i-1] + numNodes;
	for (int i=0; i<3; i++)
		for (int j=0; j<numNodes; j++)
			lengths[i][j] = 0;
			
	// and allocate a matrix that holds the state sets of the tips
	codedMatrix = new unsigned*[numTaxa];
	codedMatrix[0] = new unsigned[numTaxa * numChar];
	for (int i=1; i<numTaxa; i++)
		codedMatrix[i] = codedMatrix[i-1] + numChar;
	for (int i=0; i<numTaxa; i++)
		for (int j=0; j<numChar; j++)
			codedMatrix[i][j] = 0;

	// initialize the tip state sets
	for (int i=0; i<numTaxa; i++)
		{
		unsigned* st0 = stsPtr[0][i];
		unsigned* st1 = stsPtr[1][i];
		unsigned* st2 = stsPtr[2][i];
		for (int c=0; c<numChar; c++)
			{
            int nucCode = alignmentPtr->matrixEntry(i, c);
			int nucs[4];
			alignmentPtr->getPossibleNucs(nucCode, nucs);
			unsigned x = 0;
			for (int s=0; s<4; s++)
				{
				if ( nucs[s] == 1 )
					x += (unsigned)pow(2.0, (double)s);
				}
			codedMatrix[i][c] = x;
			(*st0) = x;
			(*st1) = x;
			(*st2) = x;
			st0++;
			st1++;
			st2++;
			}
		}
}

StateSets::~StateSets(void) {

	delete [] sts;
	delete [] stsPtr[0];
	delete [] stsPtr;
	delete [] lengths[0];
	delete [] lengths;
	delete [] codedMatrix[0];
	delete [] codedMatrix;
}

int StateSets::calcStateSetFor(int s1, int n1, int s2, int n2, int s3, int n3) {

	unsigned* ss1 = stsPtr[s1][n1];
	unsigned* ss2 = stsPtr[s2][n2];
	unsigned* ss3 = stsPtr[s3][n3];
	int numChanges = 0;
	for (int c=0; c<numChar; c++)
		{
		if ( ((*ss2) & (*ss3)) == 0 )
			{
			(*ss1) = (*ss2) | (*ss3);
			numChanges++;
			}
		else
			{
			(*ss1) = (*ss2) & (*ss3);
			}
		ss1++;
		ss2++;
		ss3++;
		}
	lengths[s1][n1] = lengths[s2][n2] + lengths[s3][n3] + numChanges;
	return numChanges;
}

int StateSets::initializeStateSets(std::vector<Node*>& dp) {

    // we assume the tree is rooted on a tip
	
	if (dp.size() == 1)
		{
		Node* p = dp[0];
		int idx = p->getIndex();
		unsigned *ssU = stsPtr[UPPER][idx];
		unsigned *ssM = stsPtr[MIDDLE][idx];
		unsigned *ssL = stsPtr[LOWER][idx];
		for (int c=0; c<numChar; c++)
			{
			unsigned x = codedMatrix[idx][c];
			(*ssU) = x;
			(*ssM) = x;
			(*ssL) = x;
			ssU++;
			ssM++;
			ssL++;
			}
		return 0;
		}

	// pass down
	int parsimonyLength = 0;
	for (int n=0; n<dp.size(); n++)
		{
		Node* p = dp[n];
		if ( p->getIsLeaf() == true && p->getAncestor() != NULL )
			{
			/* tip */
			int idx = p->getIndex();
			unsigned *ss = stsPtr[UPPER][idx];
			for (int c=0; c<numChar; c++)
				{
				(*ss) = codedMatrix[idx][c];
				ss++;
				}
			}
		else if ( p->getIsLeaf() == true && p->getAncestor() == NULL )
			{
			/* tip */
            std::vector<Node*> pDesc = p->getDescendants();
            if (pDesc.size() != 1)
                Msg::error("Expecting one descendant of tip node at root");
			int idx = p->getIndex();
			int desIdx = pDesc[0]->getIndex();
			unsigned *ss  = stsPtr[LOWER][desIdx];
			unsigned *ssU = stsPtr[UPPER][desIdx];
			for (int c=0; c<numChar; c++)
				{
				(*ss) = codedMatrix[idx][c];
				unsigned x = (*ss);
				unsigned y = (*ssU);
				unsigned zA = (x & y);
				if ( zA == 0 )
					parsimonyLength += alignmentPtr->getNumSitesOfPattern(c);
				ss++;
				ssU++;
				}
			}
		else
			{
			/* interior node */
            std::vector<Node*> pDesc = p->getDescendants();
            if (pDesc.size() != 2)
                Msg::error("Expecting two descendants of an interior node");
			unsigned* ss  = stsPtr[UPPER][p->getIndex()];
			unsigned* ssL = stsPtr[UPPER][pDesc[0]->getIndex()];
			unsigned* ssR = stsPtr[UPPER][pDesc[1]->getIndex()];
			for (int c=0; c<numChar; c++)
				{
				unsigned x = (*ssL);
				unsigned y = (*ssR);
				unsigned zA = (x & y);
				if ( zA == 0 )
					{
					unsigned zO = (x | y);
					(*ss) = zO;
					parsimonyLength += alignmentPtr->getNumSitesOfPattern(c);
					//cout << "parsimonyLength = " << parsimonyLength << endl;
					}
				else
					(*ss) = zA;
				ss++;
				ssL++;
				ssR++;
				}
			}
		}
		
	// pass up
	for (int n=(int)dp.size()-1; n>=0; n--)
		{
		Node* p = dp[n];
        Node* pAnc = p->getAncestor();
		if (pAnc != NULL)
			{
			if ( pAnc->getAncestor() != NULL )
				{
                std::vector<Node*> pAncDesc = pAnc->getDescendants();
                if (pAncDesc.size() != 2)
                    Msg::error("Expecting two descendants of an interior node");
				Node* pU;
				if (pAncDesc[0] == p)
					pU = pAncDesc[1];
				else
					pU = pAncDesc[0];
				unsigned* ss  = stsPtr[LOWER][p->getIndex()];
				unsigned* ssL = stsPtr[LOWER][pAnc->getIndex()];
				unsigned* ssU = stsPtr[UPPER][pU->getIndex()];
				for (int c=0; c<numChar; c++)
					{
					unsigned x = (*ssU);
					unsigned y = (*ssL);
					unsigned zA = (x & y);
					if ( zA == 0 )
						{
						unsigned zO = (x | y);
						(*ss) = zO;
						}
					else
						(*ss) = zA;
					ss++;
					ssU++;
					ssL++;
					}
				}
			}
		}
		
	// get state sets for center
	for (int n=0; n<dp.size(); n++)
		{
		Node* p = dp[n];
		if ( p->getAncestor() != NULL )
			{
			int idx = p->getIndex();
			unsigned* ssU = stsPtr[UPPER][idx];
			unsigned* ssL = stsPtr[LOWER][idx];
			unsigned* ssM = stsPtr[MIDDLE][idx];
			for (int c=0; c<numChar; c++)
				{
				unsigned x = (*ssU);
				unsigned y = (*ssL);
				unsigned zA = (x & y);
				if ( zA == 0 )
					{
					unsigned zO = (x | y);
					(*ssM) = zO;
					}
				else
					(*ssM) = zA;
				ssU++;
				ssL++;
				ssM++;
				}
			}
		}
		
	return parsimonyLength;
}

