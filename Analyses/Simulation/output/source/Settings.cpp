#include <iostream>
#include <cmath>
#include <string>
#include "Settings.h"



Settings::Settings(int argc, char *argv[]) {

#	if 0
	/* set up fake command-line argument string */
	char *cmdString[27];
	cmdString[ 0] = (char*)"ap";
	cmdString[ 1] = (char*)"-i";
	cmdString[ 2] = (char*)"/Users/brianmoore/Documents/Briofile/AutoParts/AutoParts/source/test.in";
	cmdString[ 3] = (char*)"-o";
	cmdString[ 4] = (char*)"/Users/brianmoore/Documents/Briofile/AutoParts/AutoParts/source/test_smart";
	cmdString[ 5] = (char*)"-l";
	cmdString[ 6] = (char*)"100";
	cmdString[ 7] = (char*)"-p";
	cmdString[ 8] = (char*)"10";
	cmdString[ 9] = (char*)"-s";
	cmdString[10] = (char*)"10";
	cmdString[11] = (char*)"-b";
	cmdString[12] = (char*)"20.0";
	cmdString[13] = (char*)"-g";
	cmdString[14] = (char*)"4";
	cmdString[15] = (char*)"-e";
	cmdString[16] = (char*)"2.0";
	cmdString[17] = (char*)"-c";
	cmdString[18] = (char*)"1";
	cmdString[19] = (char*)"-m";
	cmdString[20] = (char*)"1.1";
	cmdString[21] = (char*)"-v";
	cmdString[22] = (char*)"2.0";
	cmdString[23] = (char*)"-k";
	cmdString[24] = (char*)"3.0";
	cmdString[25] = (char*)"-d";
	cmdString[26] = (char*)"50000";
	//cmdString[23] = "-t";
	//cmdString[24] = "/Users/johnh/Working/chungking/.tree";
	argc = 27;
	argv = cmdString;
#	endif

	enum Mode { DATA_FILE, TREE_FILE, OUTPUT_FILE, BURN_IN, CHAIN_LENGTH, PRINT_FREQ, SAMPLE_FREQ, BRLEN_PARM, NUM_GAMMA_CATS, ASRV_LAMBDA, CONC_FIXED, EXP_CATS, CONC_MEAN, CONC_VAR, TUNE_T1, TUNE_T2, TUNE_T3, TUNE_T4, TUNE_T5, NONE };

	/* set default values for parameters */
	dataFilePathName       = "";
	treeFileName           = "";
	outPutFileName         = "";
	chainLength            = 1000000;
	burnIn                 = 0;
	printFrequency         = 100;
	sampleFrequency        = 100;
	brlenLambda            = 10.0;
	numGammaCats           = 4;
	asrvLambda             = 2.0;
	isConcFixed            = true;
	expNumCats             = 3.0;
	priorConcMean          = 3.0;
	priorConcVariance      = 1.0;
	tuningParm[0]          = 500.0;
	tuningParm[1]          = log(1.4);
	tuningParm[2]          = 300.0;
	tuningParm[3]          = 300.0;
	tuningParm[4]          = log(1.2);
	
	if (argc > 1)
		{
		if (argc % 2 == 0)
			{
			printUsage();
			}
			
		/* read the command-line arguments */
		int status = NONE;
		for (int i=1; i<argc; i++)
			{
			std::string cmd = argv[i];
			//std::cout << cmd << std::endl;
			if (status == NONE)
				{
				/* read the parameter specifier */
				if ( cmd == "-i" )
					status = DATA_FILE;
				else if ( cmd == "-t" )
					status = TREE_FILE;
				else if ( cmd == "-o" )
					status = OUTPUT_FILE;
				else if ( cmd == "-l" )
					status = CHAIN_LENGTH;
				else if ( cmd == "-d" )
					status = BURN_IN;
				else if ( cmd == "-p" )
					status = PRINT_FREQ;
				else if ( cmd == "-s" )
					status = SAMPLE_FREQ;
				else if ( cmd == "-b" )
					status = BRLEN_PARM;
				else if ( cmd == "-g" )
					status = NUM_GAMMA_CATS;
				else if ( cmd == "-e" )
					status = ASRV_LAMBDA;
				else if ( cmd == "-c" )
					status = CONC_FIXED;
				else if ( cmd == "-k" )
					status = EXP_CATS;
				else if ( cmd == "-m" )
					status = CONC_MEAN;
				else if ( cmd == "-v" )
					status = CONC_VAR;
				else if ( cmd == "-t1" )
					status = TUNE_T1;
				else if ( cmd == "-t2" )
					status = TUNE_T2;
				else if ( cmd == "-t3" )
					status = TUNE_T3;
				else if ( cmd == "-t4" )
					status = TUNE_T4;
				else if ( cmd == "-t5" )
					status = TUNE_T5;
				else
					{
					std::cerr << "Could not interpret option \"" << cmd << "\"." << std::endl;
					exit(1);
					}
				}
			else
				{
				/* read the parameter */
				if ( status == DATA_FILE )
					dataFilePathName = argv[i];
				else if ( status == OUTPUT_FILE )
					outPutFileName = argv[i];
				else if ( status == TREE_FILE )
					treeFileName = argv[i];
				else if ( status == CHAIN_LENGTH )
					chainLength = atoi(argv[i]);
				else if ( status == BURN_IN )
					burnIn = atoi(argv[i]);
				else if ( status == PRINT_FREQ )
					printFrequency = atoi(argv[i]);
				else if ( status == SAMPLE_FREQ )
					sampleFrequency = atoi(argv[i]);
				else if ( status == BRLEN_PARM )
					brlenLambda = atof(argv[i]);
				else if ( status == NUM_GAMMA_CATS )
					numGammaCats = atoi(argv[i]);
				else if ( status == ASRV_LAMBDA )
					asrvLambda = atof(argv[i]);
				else if ( status == CONC_FIXED )
					isConcFixed = atoi(argv[i]);
				else if ( status == EXP_CATS )
					expNumCats = atof(argv[i]);
				else if ( status == CONC_MEAN )
					priorConcMean = atof(argv[i]);
				else if ( status == CONC_VAR )
					priorConcVariance = atof(argv[i]);
				else if ( status == TUNE_T1 )
					tuningParm[0] = atof(argv[i]);
				else if ( status == TUNE_T2 )
					tuningParm[1] = atof(argv[i]);
				else if ( status == TUNE_T3 )
					tuningParm[2] = atof(argv[i]);
				else if ( status == TUNE_T4 )
					tuningParm[3] = atof(argv[i]);
				else if ( status == TUNE_T5 )
					tuningParm[4] = atof(argv[i]);
				else
					{
					std::cerr << "Unknown status reading command line information" << std::endl;
					exit(1);
					}
				status = NONE;
				}
			}
		}
	else
		{
		printUsage();
		}	

}

double Settings::getTuningParm(std::string parmNameStr) {

	if ( parmNameStr == "tree" )
		return tuningParm[0];
	else if ( parmNameStr == "asrv" )
		return tuningParm[1];
	else if ( parmNameStr == "freqs" )
		return tuningParm[2];
	else if ( parmNameStr == "rates" )
		return tuningParm[3];
	else if ( parmNameStr == "length" )
		return tuningParm[4];
	else
		{
		std::cerr << "ERROR: Unknown parameter " << parmNameStr << std::endl;
		exit(1);
		}
	return 0.0;
}

void Settings::printUsage(void) {

	std::cout << "Usage:" << std::endl;
	std::cout << "   -i : Input file name" << std::endl;
	std::cout << "   -t : Tree file name (for constraining the analysis to a fixed tree)" << std::endl;
	std::cout << "   -o : Output file name" << std::endl;
	std::cout << "   -l : Number of MCMC cycles" << std::endl;
	std::cout << "   -d : Number of MCMC cycles to discard as the \"burn-in\"" << std::endl;
	std::cout << "   -p : Print frequency" << std::endl;
	std::cout << "   -s : Sample frequency" << std::endl;
	std::cout << "   -b : Exponential parameter for branch lengths" << std::endl;
	std::cout << "   -g : Number of discrete gamma rate categories" << std::endl;
	std::cout << "   -e : Exponential parameter for shape parameter describing ASRV" << std::endl;
	std::cout << "   -c : Concentration parameter is fixed (1) or a random variable (0)" << std::endl;
	std::cout << "   -k : Prior mean of the number of categories when the concentration parameter is fixed" << std::endl;
	std::cout << "   -m : Prior mean of the concentration parameter when it is a random variable" << std::endl;
	std::cout << "   -v : Prior variance of the concentration parameter when it is a random variable" << std::endl;
	std::cout << "  -t1 : MCMC tuning parameter for the tree topology parameter" << std::endl;
	std::cout << "  -t2 : MCMC tuning parameter for the gamma shape parameter" << std::endl;
	std::cout << "  -t3 : MCMC tuning parameter for the base frequencies" << std::endl;
	std::cout << "  -t4 : MCMC tuning parameter for the substitution rates" << std::endl;
	std::cout << "  -t5 : MCMC tuning parameter for the tree length parameter" << std::endl;
	std::cout << std::endl;
	std::cout << "Example:" << std::endl;
	std::cout << "   ./AutoParts -i <input file> -o <output file>" << std::endl;
	exit(0);

}