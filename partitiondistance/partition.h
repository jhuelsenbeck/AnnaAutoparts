#ifndef PARTITION_H
#define PARTITION_H

#include <vector>

using namespace std;

class Partition {

	public:
                            Partition(int ne);
							Partition(Partition &p);
							Partition(int ne, int *p);
							Partition(vector<Partition *> &prtList);
							~Partition(void);
			    Partition   &operator=(Partition &p);
				     bool   operator==(const Partition &P);
                      int   getDegree(void) { return degree; }
					  int   getElement(int i) { return partition[i]; }
					  int   getNumElements(void) { return numElements; }
					  int   distance(Partition *p);
					 void   print(void);
	
	private:
	               double   averageDistance(vector<Partition *> &prtList);
					  int   flip(int x);
	                 void   getRgfRepresentation(void);
					 void   retreiveDiv(int whichDiv, int *div);
					 void   setDegree(void);
					 void   setElement(int i, int j);
                      int   numElements;
                      int   degree;
					  int   *partition;
	
};



#endif