#include <iostream>
#include "Alignment.h"
#include "MbRandom.h"
#include "Settings.h"
#include "Tree.h"



int main (int argc, char *argv[]) {

    // read settings
    std::string parmFile = "/Users/johnh/Desktop/autopartssim/sim_settings.ctl";
    std::string outFile  = "/Users/johnh/Desktop/autopartssim/simdata";
    Settings mySettings( parmFile );
    
    // instantiate random number object
    MbRandom myRandom;
    
    // generate a tree
    Tree myTree( mySettings.getNumTaxa(), &myRandom );
    myTree.print();
    
    // simulate the data
    Alignment mySimulation( &mySettings, &myTree, &myRandom );
    mySimulation.print();
    mySimulation.printToFile(outFile);
    
    return 0;
}
