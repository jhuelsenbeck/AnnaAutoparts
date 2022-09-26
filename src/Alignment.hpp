#ifndef Alignment_hpp
#define Alignment_hpp

#include <set>
#include <string>
#include <vector>



class Alignment {

    public:
                                    Alignment(void);
                                    Alignment(const Alignment& a);
                                    Alignment(std::string fn);
                                   ~Alignment(void);
        Alignment&                  operator+=(const Alignment& s);
        void                        compress(void);
        Alignment*                  dataForPartition(int idx);
        int                         getNumExcludedSites(void);
        int                         getNumTaxa(void) { return numTaxa; }
        std::set<int>               getPartitionIds(void) const;
        void                        getPossibleNucs(int nucCode, int* nuc);
        std::vector<std::string>    getTaxonNames(void) { return taxonNames; }
        int                         degree(void);
        int                         getNumSites(void);
        int                         getNumSitesOfPattern(int idx);
        int                         matrixEntry(int taxonIdx, int siteIdx);
        void                        print(void);
        void                        print(std::string fn);
        void                        uncompress(void);
    
    private:
        char                        convertNuc(int nucCode);
        void                        interpretString(std::string s, bool* v, int n);
        bool                        isNumber(std::string s);
        int                         longestTaxonName(void);
        int                         nucID(char nuc);
        int                         numTaxa;
        int                         numSites;
        int**                       matrix;
        bool*                       isExcluded;
        int*                        partitionId;
        std::vector<std::string>    taxonNames;
        bool                        isCompressed;
        int                         numSitePatterns;
        int*                        patternCount;
        int**                       compressedMatrix;
        int*                        compressedPartitionId;
        std::string                 pathName;
        std::string                 fileName;
        std::string                 name;
};

#endif
