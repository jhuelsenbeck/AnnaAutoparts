#ifndef Alignment_H
#define Alignment_H

#include <string>
#include <vector>

class Alignment {

	public:
                	              Alignment(std::string fileName);  
								 ~Alignment(void);
                	       void   compress(void);
                	        int   getNumTaxa(void) { return numTaxa; }
                	        int   getNumChar(void) { return (compressedData == true ? numSitePatterns : numChar); }
                	       void   getPossibleNucs (int nucCode, int *nuc);
                	        int   getNucleotide(int i, int j);
							int   getTaxonIndex(std::string ns);
							int   getPartitionId(int i) { return partitionId[i]; }
                	       void   listTaxa(void);
                	       void   print(void);
                	       void   uncompress(void);
				    std::string   getTaxonName(int i);
						    int   getNumOfPattern(int i) { return patternCount[i]; }
						   bool   getIsExcluded(int i) { return isExcluded[i]; }
						    int   getNumSubsets(void);
                	              
	private:
						   void   interpretString(std::string s, bool *v, int n);
						   bool   isNumber(std::string s);
                            int   numTaxa;
                            int   numChar;
                            int   numSitePatterns;
	   std::vector<std::string>   taxonNames;
                           bool   compressedData;
						   bool   *isExcluded;
						    int   *partitionId;
                            int   nucID(char nuc);
                            int   **matrix;
                            int   **compressedMatrix;
                            int   *patternCount;
};

#endif