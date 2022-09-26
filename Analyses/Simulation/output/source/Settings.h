#ifndef Settings_H
#define Settings_H

#include <string>


class Settings {

	public:
                            Settings(int argc, char *argv[]);
				   double   getAsrvLambda(void) { return asrvLambda; }
				   double   getBrlenLambda(void) { return brlenLambda; }
				      int   getBurnIn(void) { return burnIn; }
				      int   getChainLength(void) { return chainLength; }
			  std::string   getDataFilePathName(void) { return dataFilePathName; }
				      int   getPrintFrequency(void) { return printFrequency; }
				      int   getNumGammaCats(void) { return numGammaCats; }
			  std::string   getOutPutFileName(void) { return outPutFileName; }
					  int   getSampleFrequency(void) { return sampleFrequency; }
			  std::string   getTreeFileName(void) { return treeFileName; }
					 void   setAsrvLambda(double x) { asrvLambda = x; }
					 void   setBrlenLambda(double x) { brlenLambda = x; }
					 void   setBurnIn(int x ) { burnIn = x; }
					 void   setDataFilePathName(std::string s) { dataFilePathName = s; }
					 void   setNumGammaCats(int x) { numGammaCats = x; }
					 void   setOutPutFileName(std::string s) { outPutFileName = s; }
					 void   setTreeFileName(std::string s ) { treeFileName = s; }
					 bool   getIsConcFixed(void) { return isConcFixed; }
				   double   getExpNumCats(void) { return expNumCats; }
				   double   getPriorConcMean(void) { return priorConcMean; }
				   double   getPriorConcVariance(void) { return priorConcVariance; }
				   double   getTuningParm(std::string parmNameStr);

	private:
	                 void   printUsage(void);
				   double   brlenLambda;
				      int   chainLength;
					  int   burnIn;
			  std::string   dataFilePathName;
				      int   numGammaCats;
			  std::string   outPutFileName;
			  std::string   treeFileName;
					  int   printFrequency;
					  int   sampleFrequency;
				   double   asrvLambda;
				     bool   isConcFixed;
				   double   expNumCats;
				   double   priorConcMean;
				   double   priorConcVariance;
				   double   tuningParm[5];

};


#endif