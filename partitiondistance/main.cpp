#include <iostream>
#include <iomanip>
#include <string>
#include <istream>
#include <sstream>
#include <fstream>
#include <vector>
#include "partition.h"

using namespace std;

void readFile(string fileName, int burnIn, vector<Partition *> &partitions);
void interpretPartitionStr(string prtStr, vector<int> &prtVec);
void mySort(vector<int> &item, vector<int> &assoc, int count);
void sort2(vector<int> &item, vector<int> &assoc, int left, int right);


int main (int argc, char * const argv[]) {

//	string fileName = "/Users/brianmoore/Desktop/ap_sim1_part/ap_sim_X1_2.treelength.part";
	string fileName = "/Users/brianmoore/grunts/Grunts_k6_2.subrates.part";
	int burnIn = 2000;
	
	vector<Partition *> partitions;
	readFile(fileName, burnIn, partitions);
	
	// find the unique partitions
	vector<int> numOfThisPart;
	vector<Partition *> uniqueParts;
	vector<int> indices;
	int j = 0;
	for (vector<Partition *>::iterator p=partitions.begin(); p != partitions.end(); p++)
		{
		bool foundPart = false;
		int i = 0;
		for (vector<Partition *>::iterator q=uniqueParts.begin(); q != uniqueParts.end(); q++)
			{
			if ( (**p) == (**q) )
				{
				foundPart = true;
				break;
				}
			i++;
			}
		if ( foundPart == false )
			{
			uniqueParts.push_back( new Partition(**p) );
			numOfThisPart.push_back( 1 );
			indices.push_back( j++ );
			}
		else
			{
			numOfThisPart[i]++;
			}
		}
	mySort( numOfThisPart, indices, numOfThisPart.size() );
	
#	if 1
	int i = 0;
	for (vector<Partition *>::iterator p=partitions.begin(); p != partitions.end(); p++)
		{
		cout << setw(4) << ++i << " -- ";
		(*p)->print();
		}
		
	cout << "Unique partitions:" << endl;
	i = 0;
	double cumulative = 0.0;
	for (int i=0; i<uniqueParts.size(); i++)
		{
		double prob = (double)numOfThisPart[i] / partitions.size();
		cumulative += prob;
		cout << i + 1 << " -- " << setw(5) << numOfThisPart[i] << " " << fixed << setprecision(4) << prob << " " << cumulative << " ";
		uniqueParts[ indices[i] ]->print();
		}
	//getchar();
#	endif

	Partition avePart(partitions);
	cout << "Average partition = ";
	avePart.print();

	int *pp = new int[ partitions[0]->getNumElements() ];
	for (int i=0; i<partitions[0]->getNumElements(); i++)
		pp[i] = 0;
	int largestDegree = 0, n = 0;
	for (vector<Partition *>::iterator p=partitions.begin(); p != partitions.end(); p++)
		{
		pp[ (*p)->getDegree() ]++;
		if ( (*p)->getDegree() > largestDegree )
			largestDegree = (*p)->getDegree();
		n++;
		}
	for (int i=1; i<=largestDegree; i++)
		cout << setw(5) << i << " -- " << pp[i] << " -- " << (double)pp[i] / n << endl;
	cout << "Number of sampled partitions = " << partitions.size() << endl;
	
	/* free up memory */
	delete [] pp;
	for (vector<Partition *>::iterator p=partitions.begin(); p != partitions.end(); p++)
		delete (*p);
		
    return 0;
	
}

void readFile(string fileName, int burnIn, vector<Partition *> &partitions) {

	char lftPartDenoter = '[';
	char rhtPartDenoter = ']';
	/* open the file */
	ifstream seqStream(fileName.c_str());
	if (!seqStream) 
		{
		cerr << "Cannot open file \"" + fileName + "\"" << endl;
		exit(1);
		}

	/* read the file one line at a time */
	string linestring = "";
	string partitionStr = "";
	int line = 0;
	int numPartitions = 0;
	bool readingPartition = false;
	while( getline(seqStream, linestring).good() )
		{
		if (line > 0)
			{
			istringstream linestream(linestring);
			int ch;
			string word = "";
			do
				{
				word = "";
				linestream >> word;
				for (int i=0; i<word.size(); i++)
					{
					if ( word.at(i) == lftPartDenoter )
						{
						partitionStr = "";
						readingPartition = true;
						}
						
					if ( readingPartition == true )
						{
						if ( word.at(i) != lftPartDenoter && word.at(i) != rhtPartDenoter )
							partitionStr += word.at(i);
						}

					if ( word.at(i) == rhtPartDenoter )
						{
						readingPartition = false;
						vector<int> partitionVec;
						partitionStr += ";";
						interpretPartitionStr(partitionStr, partitionVec);
						int *tempVec = new int[partitionVec.size()];
						for (int j=0; j<partitionVec.size(); j++)
							tempVec[j] = partitionVec[j];
							
						/*for (int j=0; j<partitionVec.size(); j++)
							cout << tempVec[j] << " ";
						cout << endl;*/
						
						numPartitions++;
						if (numPartitions > burnIn)
							partitions.push_back( new Partition(partitionVec.size(), &tempVec[0]) );
						delete [] tempVec;
						}
					}

				} while ( (ch=linestream.get()) != EOF );
			}
		line++;
		}
	
	/* close the file */
	seqStream.close();

}

void interpretPartitionStr(string prtStr, vector<int> &prtVec) {
	
	string temp = "";
	for (int i=0; i<prtStr.size(); i++)
		{
		if ( prtStr.at(i) == ',' || prtStr.at(i) == ';' )
			{
			int x;
			istringstream buf(temp);
			buf >> x;
			prtVec.push_back( x );
			temp = "";
			}
		else
			temp += prtStr.at(i);
		}
}

void mySort(vector<int> &item, vector<int> &assoc, int count) {

	sort2(item, assoc, 0, count-1);
}

void sort2(vector<int> &item, vector<int> &assoc, int left, int right) {

	int i = left;
	int j = right;
	int x = item[(left+right)/2];
	do 
		{
		/*while (item[i] < x && i < right)
			i++;
		while (x < item[j] && j > left)
			j--;*/

		while (item[i] > x && i < right)
			i++;
		while (x > item[j] && j > left)
			j--;

		if (i <= j)
			{
			int y = item[i];
			item[i] = item[j];
			item[j] = y;
			
			int temp = assoc[i];
			assoc[i] = assoc[j];
			assoc[j] = temp;
			
			i++;
			j--;
			}
		} while (i <= j);
	if (left < j)
		sort2 (item, assoc, left, j);
	if (i < right)
		sort2 (item, assoc, i, right);

}

