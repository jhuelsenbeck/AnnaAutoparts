#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include "Alignment.hpp"
#include "Msg.hpp"



Alignment::Alignment(void) {

    numTaxa = 0;
    numSites = 0;
    matrix = NULL;
    isExcluded = NULL;
    partitionId = NULL;
    taxonNames.clear();
    isCompressed = false;
    numSitePatterns = 0;
    patternCount = NULL;
    compressedMatrix = NULL;
    compressedPartitionId = NULL;
    pathName = "";
    fileName = "";
    name = "Empty matrix";
}

Alignment::Alignment(const Alignment& a) {

    numTaxa = a.numTaxa;;
    numSites = a.numSites;
    
    matrix = new int*[numTaxa];
    matrix[0] = new int[numTaxa*numSites];
    for (int i=1; i<numTaxa; i++)
        matrix[i] = matrix[i-1] + numSites;
    for (int i=0; i<numTaxa; i++)
        for (int j=0; j<numTaxa; j++)
            matrix[i][j] = a.matrix[i][j];
            
    isExcluded = new bool[numSites];
    for (int i=0; i<numSites; i++)
        isExcluded[i] = a.isExcluded[i];
        
    partitionId = new int[numSites];
    for (int i=0; i<numSites; i++)
        partitionId[i] = a.partitionId[i];
        
    isCompressed = a.isCompressed;
    if (isCompressed == true)
        {
        numSitePatterns = a.numSitePatterns;
        patternCount = new int[numSitePatterns];
        for (int i=0; i<numSitePatterns; i++)
            patternCount[i] = a.patternCount[i];
        compressedMatrix = new int*[numTaxa];
        compressedMatrix[0] = new int[numTaxa*numSitePatterns];
        for (int i=1; i<numTaxa; i++)
            compressedMatrix[i] = compressedMatrix[i-1] + numSitePatterns;
        for (int i=0; i<numTaxa; i++)
            for (int j=0; j<numSitePatterns; j++)
                compressedMatrix[i][j] = a.compressedMatrix[i][j];
        compressedPartitionId = new int[numSitePatterns];
        for (int i=0; i<numSitePatterns; i++)
            compressedPartitionId[i] = a.compressedPartitionId[i];
        }
    else
        {
        numSitePatterns = 0;
        patternCount = NULL;
        compressedMatrix = NULL;
        compressedPartitionId = NULL;
        }
    
    taxonNames = a.taxonNames;
    pathName = a.pathName;
    fileName = a.fileName;
    name = a.name;
}

Alignment::Alignment(std::string fn) {

    std::cout << "   Reading data" << std::endl;
    std::cout << "   * File = \"" << fn << "\"" << std::endl << std::endl;
    
    // get last component as file name
    std::size_t found = fn.find_last_of("/\\");
    pathName = fn.substr(0,found);
    fileName = fn.substr(found+1);
    name     = fileName;

    matrix = NULL;
    numTaxa = 0;
    numSites = 0;
    isCompressed = false;
    numSitePatterns = 0;
    patternCount = NULL;
    compressedMatrix = NULL;
    compressedPartitionId = NULL;

    /* open the file */
    std::ifstream seqStream(fn.c_str());
    if (!seqStream)
        {
        std::cerr << "Cannot open file \"" + fn + "\"" << std::endl;
        exit(1);
        }

    std::string lineString = "";
    int line = 0, taxonNum = 0, pid = 0;;
    std::string theSequence = "";
    bool excludeLine = false, charSetLine = false;
    bool* tempVec = NULL;
    while( getline(seqStream, lineString).good() )
        {
        std::istringstream linestream(lineString);
        int ch;
        std::string word = "";
        int wordNum = 0;
        int siteNum = 0;
        excludeLine = false;
        charSetLine = false;
        std::string cmdString = "";
        bool foundEqualSign = false;
        do
            {
            word = "";
            linestream >> word;
            wordNum++;
            if (line == 0)
                {
                /* read the number of taxa/chars from the first line */
                int x;
                std::istringstream buf(word);
                buf >> x;
                if (wordNum == 1)
                    numTaxa = x;
                else
                    numSites = numSitePatterns = x;
                if (numTaxa > 0 && numSites > 0 && matrix == NULL)
                    {
                    matrix = new int*[numTaxa];
                    matrix[0] = new int[numTaxa * numSites];
                    for (int i=1; i<numTaxa; i++)
                        matrix[i] = matrix[i-1] + numSites;
                    for (int i=0; i<numTaxa; i++)
                        for (int j=0; j<numSites; j++)
                            matrix[i][j] = 0;
                    isExcluded = new bool[numSites];
                    partitionId = new int[numSites];
                    for (int i=0; i<numSites; i++)
                        {
                        isExcluded[i] = false;
                        partitionId[i] = -1;
                        }
                    tempVec = new bool[numSites];
                    }
                }
            else
                {
                if (wordNum == 1)
                    {
                    if ( word == "exclude" )
                        {
                        excludeLine = true;
                        foundEqualSign = false;
                        }
                    else if ( word == "charset" )
                        {
                        charSetLine = true;
                        foundEqualSign = false;
                        }
                    else
                        {
                        taxonNames.push_back(word);
                        taxonNum++;
                        }
                    }
                else
                    {
                    if (excludeLine == true || charSetLine == true)
                        {
                        for (int i=0; i<word.size(); i++)
                            {
                            if (word.at(i) == '=')
                                foundEqualSign = true;
                            if (foundEqualSign == true)
                                {
                                if ( isdigit(word.at(i)) )
                                    cmdString += word.at(i);
                                else if ( word.at(i) == '-' )
                                    cmdString += " - ";
                                else if (word.at(i) == '\\')
                                    cmdString += " \\ ";
                                }
                            }
                        cmdString += " ";
                        for (int i=0; i<word.size(); i++)
                            {
                            if ( word.at(i) == ';' )
                                {
                                interpretString(cmdString, tempVec, numSites);
                                if (excludeLine == true)
                                    {
                                    for (int i=0; i<numSites; i++)
                                        if (tempVec[i] == true)
                                            isExcluded[i] = true;
                                    }
                                else if (charSetLine == true)
                                    {
                                    pid++;
                                    for (int i=0; i<numSites; i++)
                                        if (tempVec[i] == true)
                                            partitionId[i] = pid;
                                    }
                                }
                            }
                        }
                    else
                        {
                        for (int i=0; i<word.length(); i++)
                            {
                            char site = word.at(i);
                            matrix[taxonNum-1][siteNum++] = nucID(site);
                            }
                        }
                    }
                }
            } while ( (ch=linestream.get()) != EOF );
            
        // NOTE: We probably do not need this bit of code.
        if (line == 0)
            {
            /* the first line should contain the number of taxa and the sequence length */
            std::istringstream buf(word);
            //buf >> genomeSize;
            }
        else
            {
            for (int i=0; i<word.length(); i++)
                {
                char site = word.at(i);
                if (tolower(site) == 'a' || tolower(site) == 'c' || tolower(site) == 'g' || tolower(site) == 't')
                    theSequence += tolower(site);
                }
            }
        line++;
        }

    delete [] tempVec;
    
    /* close the file */
    seqStream.close();

}

Alignment::~Alignment(void) {

    delete [] matrix[0];
    delete [] matrix;
    delete [] isExcluded;
    delete [] partitionId;
    uncompress();
}

Alignment& Alignment::operator+=(const Alignment& rhs) {

    std::cout << "concatenating " << this->name << " with " << rhs.name << std::endl;
    
    // check that the two matrices can be concatenated
    if (this->numTaxa == 0)
        {
        this->numTaxa = rhs.numTaxa;
        this->taxonNames = rhs.taxonNames;
        }
    
    std::cout << "This numTaxa = " << this->taxonNames.size() << std::endl;
    std::cout << "RHS numTaxa = " << rhs.taxonNames.size() << std::endl;
    for (int i=0; i<this->taxonNames.size(); i++)
        std::cout << i << " " << this->taxonNames[i] << std::endl;
    for (int i=0; i<rhs.taxonNames.size(); i++)
        std::cout << i << " " << rhs.taxonNames[i] << std::endl;
        
    if (this->numTaxa != rhs.numTaxa)
        Msg::error("The matrices have different taxon membership and cannot be concatenated");
    for (int i=0; i<numTaxa; i++)
        {
        if (this->taxonNames[i] != rhs.taxonNames[i])
            Msg::error("The matrices have different taxa and cannot be concatenated");
        }
        
    // uncompress this data if it is compressed
    bool dataWasOriginallyCompressed = false;
    if (this->isCompressed == true)
        {
        this->uncompress();
        dataWasOriginallyCompressed = true;
        }
    
    // reallocate the matrix
    int** tempMatrix = new int*[numTaxa];
    tempMatrix[0] = new int[numTaxa * (this->numSites + rhs.numSites)];
    for (int i=1; i<numTaxa; i++)
        tempMatrix[i] = tempMatrix[i-1] + (this->numSites + rhs.numSites);
    for (int i=0; i<numTaxa; i++)
        for (int j=0; j<(this->numSites + rhs.numSites); j++)
            tempMatrix[i][j] = 0;
            
    // fill it out
    for (int i=0; i<numTaxa; i++)
        {
        int k = 0;
        for (int j=0; j<this->numSites; j++)
            tempMatrix[i][k++] = this->matrix[i][j];
        for (int j=0; j<rhs.numSites; j++)
            tempMatrix[i][k++] = rhs.matrix[i][j];
        }
    if (this->matrix != NULL)
        {
        delete [] this->matrix[0];
        delete [] this->matrix;
        }
    this->matrix = tempMatrix;
    this->matrix[0] = tempMatrix[0];
    
    // fill out the partiton information
    std::set<int> lhsParts = this->getPartitionIds();
    std::set<int> rhsParts = rhs.getPartitionIds();
    std::map<int,int> lhsMap;
    std::map<int,int> rhsMap;
    int k = 1;
    for (int key : lhsParts)
        {
        lhsMap.insert( std::make_pair(key,k) );
        k++;
        }
    for (int key : rhsParts)
        {
        rhsMap.insert( std::make_pair(key,k) );
        k++;
        }
    
    int* tempPartitionId = new int[this->numSites + rhs.numSites];
    for (int i=0; i<this->numSites; i++)
        {
        int key = this->partitionId[i];
        std::map<int,int>::iterator it = lhsMap.find(key);
        if (it == lhsMap.end())
            Msg::error("Could not find partition id in map");
        tempPartitionId[i] = it->second;
        }
    for (int i=this->numSites, k=0; i<this->numSites+rhs.numSites; i++)
        {
        int key = rhs.partitionId[k++];
        std::map<int,int>::iterator it = rhsMap.find(key);
        if (it == rhsMap.end())
            Msg::error("Could not find partition id in map");
        tempPartitionId[i] = it->second;
        }
    if (this->partitionId != NULL)
        delete [] this->partitionId;
    this->partitionId = tempPartitionId;
    
    // fill out the isExcluded vector
    bool* tempIsExcluded = new bool[this->numSites + rhs.numSites];
    for (int i=0; i<this->numSites; i++)
        tempIsExcluded[i] = this->isExcluded[i];
    for (int i=this->numSites, k=0; i<this->numSites+rhs.numSites; i++)
        tempIsExcluded[i] = rhs.isExcluded[k++];
    if (this->isExcluded != NULL)
        delete [] this->isExcluded;
    this->isExcluded = tempIsExcluded;
            
    // finally, change the number of sites
    this->numSites += rhs.numSites;
    
    // deal with compresseion
    if (dataWasOriginallyCompressed == true)
        this->compress();

    return *this;
}

void Alignment::compress(void) {

    if (isCompressed == true)
        return;
    
    // find the unique site patterns and count them up
    int* tempPatternIds = new int[numSites];
    for (int i=0; i<numSites; i++)
        tempPatternIds[i] = -1;
    numSitePatterns = 0;
    for (int i=0; i<numSites; i++)
        {
        if (tempPatternIds[i] == -1 && isExcluded[i] == false)
            {
            tempPatternIds[i] = numSitePatterns;
            for (int j=i+1; j<numSites; j++)
                {
                if (tempPatternIds[j] == -1)
                    {
                    bool isSame = true;
                    for (int k=0; k<numTaxa; k++)
                        {
                        if (matrix[k][i] != matrix[k][j])
                            {
                            isSame = false;
                            break;
                            }
                        }
                    // only equate the two patterns if the sites are the same, they have the same
                    // partition id, and they are not excluded
                    if (isSame == true && partitionId[i] == partitionId[j] && isExcluded[j] == false)
                        tempPatternIds[j] = numSitePatterns;
                    }
                }
            numSitePatterns++;
            }
        }
#   if 0
    for (int j=0; j<numSites; j++)
        {
        std::cout << std::setw(5) << j << " -- " << std::setw(5) << tempPatternIds[j] << " -- ";
        for (int i=0; i<numTaxa; i++)
            std::cout << convertNuc(matrix[i][j]) << " ";
        std::cout << std::endl;
        }
#   endif

    // the number of sites that have not been compressed should be equal to the number of excluded sites. Check!
    int numNotAccountedFor = 0;
    for (int i=0; i<numSites; i++)
        {
        if (tempPatternIds[i] == -1)
            numNotAccountedFor++;
        }
    if (numNotAccountedFor != getNumExcludedSites())
        Msg::error("Expecting the number of sites not compressed to be equal to the number of excluded sites");
        
    // fill in the compressed data matrix
    if (compressedMatrix != NULL)
        Msg::error("Expecting compressedMatrix to be NULL");
    if (patternCount != NULL)
        Msg::error("Expecting patternCount to be NULL");
    if (compressedPartitionId != NULL)
        Msg::error("Expecting compressedPartitionId to be NULL");
    compressedMatrix = new int*[numTaxa];
    compressedMatrix[0] = new int[numTaxa*numSitePatterns];
    for (int i=1; i<numTaxa; i++)
        compressedMatrix[i] = compressedMatrix[i-1] + numSitePatterns;
    patternCount = new int[numSitePatterns];
    compressedPartitionId = new int[numSitePatterns];
    
    for (int j=0; j<numSitePatterns; j++)
        {
        int numOfThisPattern = 0;
        bool foundPattern = false;
        for (int k=0; k<numSites; k++)
            {
            if (tempPatternIds[k] == j)
                {
                numOfThisPattern++;
                if (foundPattern == false)
                    {
                    for (int n=0; n<numTaxa; n++)
                        compressedMatrix[n][j] = matrix[n][k];
                    }
                foundPattern = true;
                compressedPartitionId[j] = partitionId[k];
                }
            }
        if (foundPattern == false)
            Msg::error("Could not find pattern");
        patternCount[j] = numOfThisPattern;
        }
        
    // check that the number of sites is equal to the number of non-excluded sites
    int n = 0;
    for (int i=0; i<numSitePatterns; i++)
        n += patternCount[i];
    if (n != numSites - getNumExcludedSites())
        Msg::error("Problem compressing data");
    
    isCompressed = true;
    
    delete [] tempPatternIds;
}

char Alignment::convertNuc(int nucCode) {

    if (nucCode == 1)
        return 'A';
    else if (nucCode == 2)
        return 'C';
    else if (nucCode == 3)
        return 'M';
    else if (nucCode == 4)
        return 'G';
    else if (nucCode == 5)
        return 'R';
    else if (nucCode == 6)
        return 'S';
    else if (nucCode == 7)
        return 'V';
    else if (nucCode == 8)
        return 'T';
    else if (nucCode == 9)
        return 'W';
    else if (nucCode == 10)
        return 'Y';
    else if (nucCode == 11)
        return 'H';
    else if (nucCode == 12)
        return 'K';
    else if (nucCode == 13)
        return 'D';
    else if (nucCode == 14)
        return 'B';
    else if (nucCode == 15)
        return 'N';
    else if (nucCode == 16)
        return '-';
    else
        Msg::error("Unknown character code");
    return ' ';
}

Alignment* Alignment::dataForPartition(int idx) {
        
    // check that partition id exists
    bool idExists = false;
    for (int i=0; i<numSites; i++)
        {
        if (partitionId[i] == idx)
            {
            idExists = true;
            break;
            }
        }
    if (idExists == false)
        Msg::error("No sites of this pattern index");
    
    // count the number of sites with this id
    int n = 0;
    for (int i=0; i<numSites; i++)
        {
        if (partitionId[i] == idx)
            n++;
        }
            
    // make a new alignment object
    Alignment* aln = new Alignment;
    aln->numSites = n;
    aln->numTaxa = this->numTaxa;
    aln->taxonNames = this->taxonNames;
    
    // set its name
    aln->pathName = this->pathName;
    aln->fileName = this->fileName;
    aln->name     = this->name + " (Partition ID " + std::to_string(idx) + ")";

    aln->matrix = new int*[aln->numTaxa];
    aln->matrix[0] = new int[aln->numTaxa * aln->numSites];
    for (int i=1; i<aln->numTaxa; i++)
        aln->matrix[i] = aln->matrix[i-1] + aln->numSites;
    for (int i=0; i<aln->numTaxa; i++)
        for (int j=0; j<aln->numSites; j++)
            aln->matrix[i][j] = 0;
    aln->isExcluded = new bool[aln->numSites];
    aln->partitionId = new int[aln->numSites];
    
    // fill in the observations for this new object
    n = 0;
    for (int i=0; i<this->numSites; i++)
        {
        if (this->partitionId[i] == idx)
            {
            for (int j=0; j<numTaxa; j++)
                aln->matrix[j][n] = this->matrix[j][i];
            aln->isExcluded[n] = this->isExcluded[i];
            aln->partitionId[n] = idx;
            n++;
            }
        }

    aln->isCompressed = false;
    aln->patternCount = NULL;
    aln->compressedMatrix = NULL;
    aln->compressedPartitionId = NULL;
    if (this->isCompressed == true)
        aln->compress();

    return aln;
}

int Alignment::degree(void) {

    std::set<int> pids;
    for (int i=0; i<numSites; i++)
        pids.insert(partitionId[i]);
    return (int)pids.size();
}

int Alignment::getNumSites(void) {

    if (isCompressed == false)
        return numSites;
    else
        return numSitePatterns;
}

int Alignment::getNumSitesOfPattern(int idx) {

    if (isCompressed == false)
        return 1;
    else
        {
        if (idx < numSitePatterns)
            return patternCount[idx];
        Msg::error("Site index is too large");
        }
    return 0;
}

std::set<int> Alignment::getPartitionIds(void) const {

    std::set<int> pids;
    for (int i=0; i<numSites; i++)
        pids.insert(partitionId[i]);
    return pids;
}

int Alignment::getNumExcludedSites(void) {

    int n = 0;
    for (int i=0; i<numSites; i++)
        {
        if (isExcluded[i] == true)
            n++;
        }
    return n;
}

void Alignment::interpretString(std::string s, bool* v, int n) {

    for (int i=0; i<n; i++)
        v[i] = false;
    int nums[3];
    int numToRemember = 0;
    //(*log) << "s = \"" << s << "\"" << std::endl;

    /* push the individual words (numbers, hyphens, or back slashes into a vector */
    std::vector<std::string> cmds;
    std::istringstream linestream(s);
    int ch;
    std::string word = "";
    do
        {
        word = "";
        linestream >> word;
        if (word != "")
            {
            cmds.push_back( word );
            }
        } while( (ch=linestream.get()) != EOF );
    
    /* loop over the vector of individual words and correctly interpret things */
    int i = 0;
    for (std::vector<std::string>::iterator p=cmds.begin(); p != cmds.end(); p++)
        {
        //(*log) << "\"" << (*p) << "\"" << std::endl;
        if ( isdigit((*p)[0]) )
            {
            
            /* we can potentially complete a sentence */
            std::string prevWord = "";
            if (i > 0)
                prevWord = cmds[i-1];
            std::string nextWord = "";
            if (i != cmds.size() - 1)
                nextWord = cmds[i+1];
            
            int x;
            std::istringstream buf(cmds[i]);
            buf >> x;

            if (prevWord == "" || isNumber(prevWord) == true)
                {
                nums[0] = x;
                numToRemember = 1;
                }
            else if (prevWord == "-")
                {
                nums[1] = x;
                numToRemember = 2;
                }
            else if (prevWord == "\\")
                {
                nums[2] = x;
                numToRemember = 3;
                }
            else
                {
                std::cerr << "ERROR: Problem interpreting string" << std::endl;
                exit(1);
                }
            
            if ( (prevWord == "" || isNumber(prevWord) == true) && (nextWord == "" || isNumber(nextWord) == true) )
                {
                v[nums[0]-1] = true;
                numToRemember = 0;
                }
            else if ( prevWord == "-" && (nextWord == "" || isNumber(nextWord) == true) )
                {
                for (int i=nums[0]-1; i<nums[1]; i++)
                    v[i] = true;
                numToRemember = 0;
                }
            else if ( prevWord == "\\" )
                {
                for (int i=nums[0]-1, k=nums[2]; i<nums[1]; i++, k++)
                    if ( k % nums[2] == 0 )
                        v[i] = true;
                numToRemember = 0;
                }
            }
        i++;
        }
}

bool Alignment::isNumber(std::string s) {

    if (s == "")
        return false;

    bool isnum = true;
    for (int i=0; i<s.size(); i++)
        if ( !isdigit(s[i]) )
            isnum = false;
    return isnum;
}

int Alignment::longestTaxonName(void) {

    int len = 0;
    for (int i=0; i<taxonNames.size(); i++)
        {
        if (taxonNames[i].length() > len)
            len = (int)taxonNames[i].length();
        }
    return len;
}

int Alignment::matrixEntry(int taxonIdx, int siteIdx) {

    if (isCompressed == false)
        {
        if (taxonIdx >= numTaxa || siteIdx >= numSites)
            Msg::error("Matrix index out of bounds");
        return matrix[taxonIdx][siteIdx];
        }
    else
        {
        if (taxonIdx >= numTaxa || siteIdx >= numSitePatterns)
            Msg::error("Compressed matrix index out of bounds");
        return compressedMatrix[taxonIdx][siteIdx];
        }
}

/*-------------------------------------------------------------------
|
|   NucID:
|
|   Take a character, nuc, and return an integer:
|
|       nuc        returns
|        A            1
|        C            2
|        G            4
|        T/U          8
|        R            5
|        Y           10
|        M            3
|        K           12
|        S            6
|        W            9
|        H           11
|        B           14
|        V            7
|        D           13
|        N/?         15
|        -           16
|
-------------------------------------------------------------------*/
int Alignment::nucID(char nuc) {

	char		n;
	
	if (nuc == 'U' || nuc == 'u')
		n = 'T';
	else
		n = nuc;

	if (n == 'A' || n == 'a')
		{
		return 1;
		}
	else if (n == 'C' || n == 'c')
		{
		return 2;
		}
	else if (n == 'G' || n == 'g')
		{
		return 4;
		}
	else if (n == 'T' || n == 't')
		{
		return 8;
		}
	else if (n == 'R' || n == 'r')
		{
		return 5;
		}
	else if (n == 'Y' || n == 'y')
		{
		return 10;
		}
	else if (n == 'M' || n == 'm')
		{
		return 3;
		}
	else if (n == 'K' || n == 'k')
		{
		return 12;
		}
	else if (n == 'S' || n == 's')
		{
		return 6;
		}
	else if (n == 'W' || n == 'w')
		{
		return 9;
		}
	else if (n == 'H' || n == 'h')
		{
		return 11;
		}
	else if (n == 'B' || n == 'b')
		{
		return 14;
		}
	else if (n == 'V' || n == 'v')
		{
		return 7;
		}
	else if (n == 'D' || n == 'd')
		{
		return 13;
		}
	else if (n == 'N' || n == 'n')
		{
		return 15;
		}
	else if (n == '?')
		{
		return 15;
		}
	else if (n == '-')
		{
		return 16;
		}
	else
		return -1;
}

/*-------------------------------------------------------------------
|
|   GetPossibleNucs:
|
|   This function initializes a vector, nuc[MAX_NUM_STATES]. The four elements
|   of nuc correspond to the four nucleotides in alphabetical order.
|   We are assuming that the nucCode is a binary representation of
|   the nucleotides that are consistent with the observation. For
|   example, if we observe an A, then the nucCode is 1 and the
|   function initalizes nuc[0] = 1 and the other elements of nuc
|   to be 0.
|
|   Observation    nucCode        nuc
|        A            1           1000
|        C            2           0100
|        G            4           0010
|        T            8           0001
|        R            5           1010
|        Y           10           0101
|        M            3           1100
|        K           12           0011
|        S            6           0110
|        W            9           1001
|        H           11           1101
|        B           14           0111
|        V            7           1110
|        D           13           1011
|        N - ?       15           1111
|
-------------------------------------------------------------------*/
void Alignment::getPossibleNucs(int nucCode, int* nuc) {

	if (nucCode == 1)
		{
		nuc[0] = 1;
		nuc[1] = 0;
		nuc[2] = 0;
		nuc[3] = 0;
		}
	else if (nucCode == 2)
		{
		nuc[0] = 0;
		nuc[1] = 1;
		nuc[2] = 0;
		nuc[3] = 0;
		}
	else if (nucCode == 3)
		{
		nuc[0] = 1;
		nuc[1] = 1;
		nuc[2] = 0;
		nuc[3] = 0;
		}
	else if (nucCode == 4)
		{
		nuc[0] = 0;
		nuc[1] = 0;
		nuc[2] = 1;
		nuc[3] = 0;
		}
	else if (nucCode == 5)
		{
		nuc[0] = 1;
		nuc[1] = 0;
		nuc[2] = 1;
		nuc[3] = 0;
		}
	else if (nucCode == 6)
		{
		nuc[0] = 0;
		nuc[1] = 1;
		nuc[2] = 1;
		nuc[3] = 0;
		}
	else if (nucCode == 7)
		{
		nuc[0] = 1;
		nuc[1] = 1;
		nuc[2] = 1;
		nuc[3] = 0;
		}
	else if (nucCode == 8)
		{
		nuc[0] = 0;
		nuc[1] = 0;
		nuc[2] = 0;
		nuc[3] = 1;
		}
	else if (nucCode == 9)
		{
		nuc[0] = 1;
		nuc[1] = 0;
		nuc[2] = 0;
		nuc[3] = 1;
		}
	else if (nucCode == 10)
		{
		nuc[0] = 0;
		nuc[1] = 1;
		nuc[2] = 0;
		nuc[3] = 1;
		}
	else if (nucCode == 11)
		{
		nuc[0] = 1;
		nuc[1] = 1;
		nuc[2] = 0;
		nuc[3] = 1;
		}
	else if (nucCode == 12)
		{
		nuc[0] = 0;
		nuc[1] = 0;
		nuc[2] = 1;
		nuc[3] = 1;
		}
	else if (nucCode == 13)
		{
		nuc[0] = 1;
		nuc[1] = 0;
		nuc[2] = 1;
		nuc[3] = 1;
		}
	else if (nucCode == 14)
		{
		nuc[0] = 0;
		nuc[1] = 1;
		nuc[2] = 1;
		nuc[3] = 1;
		}
	else if (nucCode == 15)
		{
		nuc[0] = 1;
		nuc[1] = 1;
		nuc[2] = 1;
		nuc[3] = 1;
		}
	else if (nucCode == 16)
		{
		nuc[0] = 1;
		nuc[1] = 1;
		nuc[2] = 1;
		nuc[3] = 1;
		}
}

void Alignment::print(void) {

    std::cout << "Data for alignment " << name << std::endl;
    std::cout << "*  File name               = " << fileName << std::endl;
    std::cout << "*  Path to file            = " << pathName << std::endl;
    std::cout << "*  Number of taxa          = " << numTaxa << std::endl;
    if (isCompressed == false)
        {
        std::cout << "*  Number of sites         = " << getNumSites() << std::endl;
        }
    else
        {
        int ns = 0;
        for (int c=0; c<numSitePatterns; c++)
            ns += patternCount[c];
        std::cout << "*  Number of sites         = " << ns << std::endl;
        std::cout << "*  Number of site patterns = " << numSitePatterns << std::endl;
        }
    std::cout << "*  Number of partitions    = " << getPartitionIds().size() << std::endl;
    if (isCompressed == true)
        std::cout << "*  Is data compressed      = true" << std::endl;
    else
        std::cout << "*  Is data compressed      = false" << std::endl;

    if (isCompressed == false)
        {
        for (int c=0; c<numSites; c++)
            {
            std::cout << std::setw(5) << c << " -- ";
            std::cout << std::setw(2) << partitionId[c] << " -- ";
            for (int t=0; t<numTaxa; t++)
                {
                std::cout << std::setw(2) << convertNuc(matrix[t][c]) << " ";
                }
            std::cout << std::endl;
            }
        }
    else
        {
        for (int c=0; c<numSitePatterns; c++)
            {
            std::cout << std::setw(5) << c << " -- ";
            std::cout << std::setw(2) << compressedPartitionId[c] << " -- ";
            std::cout << std::setw(5) << patternCount[c] << " -- ";
            for (int t=0; t<numTaxa; t++)
                {
                std::cout << std::setw(2) << convertNuc(compressedMatrix[t][c]) << " ";
                }
            std::cout << std::endl;
            }
        }
}

void Alignment::print(std::string fn) {

    std::ofstream nexOut;
    nexOut.open( fn.c_str(), std::ios::out );
    
    nexOut << "#NEXUS\n\n";
    nexOut << "begin data;\n";
    nexOut << "   dimensions ntax=" << numTaxa << " nchar=" << numSites << ";\n";
    nexOut << "   format datatype=dna gap=-;\n";
    nexOut << "   matrix\n";
    int len = longestTaxonName();
    for (int i=0; i<numTaxa; i++)
        {
        nexOut << "   " << taxonNames[i];
        for (int j=0; j<len-taxonNames[i].length(); j++)
            nexOut << " ";
        nexOut << "  ";
        for (int j=0; j<numSites; j++)
            {
            nexOut << convertNuc(matrix[i][j]);
            }
        }
    nexOut << "   ;\n";
    nexOut << "end;\n";
    
    nexOut.close();
}

void Alignment::uncompress(void) {

    if (isCompressed == false)
        return;
        
    numSitePatterns = 0;
    delete [] patternCount;
    patternCount = NULL;
    delete [] compressedMatrix[0];
    delete [] compressedMatrix;
    compressedMatrix = NULL;
    delete [] compressedPartitionId;
    compressedPartitionId = NULL;
}
