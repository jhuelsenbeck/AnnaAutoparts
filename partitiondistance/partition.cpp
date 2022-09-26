#include "partition.h"
#include "hungarian.h"
#include "MbRandom.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <iterator>

using namespace std;



Partition::Partition(int ne) {

	numElements = ne;
	partition =  new int[numElements];
	for (int i=0; i<numElements; i++)
		partition[i] = 1;
	getRgfRepresentation();
	setDegree();

}

Partition::Partition(Partition &p) {

	numElements = p.numElements;
	partition =  new int[numElements];
	for (int i=0; i<numElements; i++)
		partition[i] = p.partition[i];
	getRgfRepresentation();
	setDegree();

}

Partition::Partition(int ne, int *p) {

	numElements = ne;
	partition =  new int[numElements];
	for (int i=0; i<numElements; i++)
		partition[i] = p[i];
	getRgfRepresentation();
	setDegree();
	
}

Partition::Partition(vector<Partition *> &prtList) {

	/* check that all of the partitions in the list have the
	   same number of elements */
	numElements = prtList[0]->getNumElements();
	for (vector<Partition *>::iterator p=prtList.begin(); p != prtList.end(); p++)
		{
		int x = (*p)->getNumElements();
		if (x != numElements)
			{
			cerr << "ERROR: The partitions differ in size" << endl;
			exit(1);
			}
		}
	
	/* heuristic search looking for better partitions (i.e., partitions that
	   minimize the squared distance to all of the partitions in the
	   list of partitions, prtList) */
	MbRandom myRandom;
	for (int i=0; i<1000; i++)
		myRandom.uniformRv();
	int whichPartition = (int)(myRandom.uniformRv() * prtList.size());
	Partition *avePart[2];
	avePart[0] = new Partition( *prtList[whichPartition] );
	avePart[1] = new Partition( *prtList[whichPartition] );
	int curAve = 0;
	
	double bestDist = avePart[0]->averageDistance(prtList);
	avePart[0]->print();
	avePart[1]->print();
	
	bool foundBetter = false;
	do
		{
		foundBetter = false;
		for (int i=0; i<numElements; i++)
			{
			int newAve = flip(curAve);
			int deg = avePart[newAve]->getDegree();
			for (int j=1; j<=deg+1; j++)
				{
				avePart[newAve]->setElement(i, j);
				avePart[newAve]->getRgfRepresentation();
				avePart[newAve]->setDegree();
				double d = avePart[newAve]->averageDistance(prtList);
#				if 1
				cout << setw(5) << i << setw(5) << numElements << setw(5) << j << " " << fixed << setprecision(4) << d << "  " << bestDist << " -- ";
				avePart[newAve]->print();
#				endif
				if ( d < bestDist )
					{
					curAve = newAve;
					newAve = flip(curAve);
					bestDist = d;
					foundBetter = true;
#					if 0
					cout << setw(5) << avePart[curAve]->getDegree() << " " << fixed << setprecision(4) << bestDist << " -- ";
					avePart[curAve]->print();
#					endif
					}
				*avePart[newAve] = *avePart[curAve];
				
				}
			}
		} while (foundBetter == true);
	
	//cout << "Average partition (" << avePart[curAve]->averageDistance(prtList) << ") = ";
	//avePart[curAve]->print();

	/* initialize the partition with the contents of the average partition */
	partition =  new int[numElements];
	for (int i=0; i<numElements; i++)
		partition[i] = avePart[curAve]->getElement( i );
	getRgfRepresentation();
	setDegree();

}

Partition::~Partition(void) {

	delete [] partition;
}

Partition &Partition::operator=(Partition &p) {

	if (this != &p)
		{
		if (numElements != p.numElements)
			{
			delete [] partition;
			partition = new int[p.numElements];
			}
		numElements = p.numElements;
		degree = p.degree;
		for (int i=0; i<numElements; i++)
			partition[i] = p.partition[i];
		}
	return *this;
}

bool Partition::operator==(const Partition &p) {

	if (numElements != p.numElements)
		return false;
	if (degree != p.degree)
		return false;
	for (int i=0; i<numElements; i++)
		if (partition[i] != p.partition[i])
			return false;
	return true;
}

void Partition::print(void) {

	cout << "(";
	for (int i=0; i<numElements; i++)
		{
		cout << partition[i];
		if ( i + 1 < numElements )
			cout << ",";
		}
	cout << ") [" << numElements << "," << degree << "]" << endl;
}

void Partition::setDegree(void) {

	/* Find the largest number in the partition. This number will be the
	   degree of the model. This calculation assumes that the partition
	   is in the restricted growth function (RGF) format. */

	int largestNumFound = 0;
	int smallestNumFound = partition[0];
	for (int i=0; i<numElements; i++)
		{
		if (partition[i] > largestNumFound)
			largestNumFound = partition[i];
		}
	degree = largestNumFound - smallestNumFound + 1;
}

int Partition::distance(Partition *p) {

	Partition *p1 = this;
	Partition *p2 = p;
	
	if ( p1->getNumElements() != p2->getNumElements() )
		return -1;
		
	int *prt1 = &p1->partition[0];
	int *prt2 = &p2->partition[0];
	
	int *div1 = new int[numElements];
	int *div2 = new int[numElements];
	
	int d1 = p1->getDegree();
	int d2 = p2->getDegree();

	int nrows, ncols;
	int *t;
	if (d1 > d2)
		{
		ncols = d1;
		nrows = d2;
		}
	else
		{
		ncols = d2;
		nrows = d1;
		t = prt1;
		prt1 = prt2;
		prt2 = t;
		p1 = p;
		p2 = this;
		}

	int **r = new int*[nrows];
	r[0] = new int[nrows * ncols];
	for (int i=1; i<nrows; i++)
		r[i] = r[i-1] + ncols;
	for (int i=0; i<nrows; i++)
		for (int j=0; j<ncols; j++)
			r[i][j] = 0;
	for (int i=0; i<nrows; i++)
		{
		p2->retreiveDiv(i+1, div2);
		for (int j=0; j<ncols; j++)
			{
			p1->retreiveDiv(j+1, div1);
			for (int k=0; k<numElements; k++)
				if (div1[k] == 1 && div2[k] == 1)
					r[i][j]++;
			}
		}
#	if 0
	for (int i=0; i<nrows; i++)
		{
		for (int j=0; j<ncols; j++)
			std::cout << r[i][j] << ",";
		std::cout << std::endl;
		}
	std::cout << "done" << std::endl;
	getchar();
#	endif
	Hungarian *h = new Hungarian(r, nrows, ncols);
	int d = numElements - h->getAssignmentCost();
	

	delete [] r[0];
	delete [] r;
	delete [] div1;
	delete [] div2;

	return d;
	
}

double Partition::averageDistance(vector<Partition *> &prtList) {

	double sum = 0.0;
	int n = 0;
	for (vector<Partition *>::iterator p=prtList.begin(); p != prtList.end(); p++)
		{
		int d = distance( *p );
		sum += (double)(d * d);
		n++;
		}
	return sum / n;
	
}

int Partition::flip(int x) {

	if (x == 0)
		return 1;
	return 0;
	
}

void Partition::getRgfRepresentation(void) {

	bool *x = new bool[numElements];
	int *assignment = new int[numElements];
	for (int i=0; i<numElements; i++)
		{
		x[i] = false;
		assignment[i] = partition[i];
		partition[i] = -1;
		}
		
	for (int i=0, idx=1; i<numElements; i++)
		{
		if (x[i] == false)
			{
			for (int j=0; j<numElements; j++)
				{
				if (x[j] == false && assignment[j] == assignment[i])
					{
					partition[j] = idx;
					x[j] = true;
					}
				}
			idx++;
			}
		}
	
	delete [] x;
	delete [] assignment;

}

void Partition::retreiveDiv(int whichDiv, int *div) {

	for (int i=0; i<numElements; i++)
		div[i] = 0;
		
	for (int i=0; i<numElements; i++)
		{
		if (partition[i] == whichDiv)
			div[i] = 1;
		}

}

void Partition::setElement(int i, int j) {

	partition[i] = j;
	
}