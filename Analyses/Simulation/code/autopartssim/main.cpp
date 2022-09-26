#include <iostream>
#include <vector>
#include "Alignment.h"
#include "MbRandom.h"
#include "Settings.h"
#include "Tree.h"



int main (int argc, char *argv[]) {

    // read settings
    std::string parmFile = "/Users/brianmoore/Documents/Briofile/AutoParts/Analyses/Simulation/autopartssim/sim_settings.ctl";
    std::string outFile  = "/Users/brianmoore/Documents/Briofile/AutoParts/Analyses/Simulation/autopartssim/sim";
    Settings mySettings( parmFile );
    
    // instantiate random number object
    MbRandom myRandom;
    
    // generate a tree
    Tree myTree( mySettings.getNumTaxa(), &myRandom );
    myTree.print();
    
    // simulate the data
	std::vector<Alignment*> simulations;
	for (int i=0; i<mySettings.getNumReplicates(); i++)
		{
		std::cout << "Simulating replicate number " << i+1 << std::endl;
		Alignment* a = new Alignment( &mySettings, &myTree, &myRandom );
		simulations.push_back( a );
		}
	for (int i=0; i<mySettings.getNumReplicates(); i++)
		{
		//simulations[i]->print();
		simulations[i]->printToFile(outFile, i);
		}
		
	for (int i=0; i<mySettings.getNumReplicates(); i++)
		delete simulations[i];
		
    return 0;
}
