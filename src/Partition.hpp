#ifndef Partition_H
#define Partition_H

#include <map>
#include <set>
#include <string>
#include <vector>
class Subset;


class Partition {

    friend class PartitionLessThan;
    
    public:
                                Partition(std::vector<int> x);
                                Partition(const Partition& p);
                                Partition(std::map<Partition*,int>& partList, double& score);
                               ~Partition(void);
        Partition&              operator=(Partition& p);
        int                     degree(void) { return (int)subsets.size(); }
        int                     degree(void) const { return (int)subsets.size(); }
        int                     distance(Partition* p);
        Subset*                 findSubsetWithElement(int x);
        std::vector<bool>&      getBitRepresentationForSubset(int idx) const;
        int                     getNumElements(void) { return numElements; }
        std::string             getRgfString(void);
        std::string             getRgfString(void) const;
        void                    print(void);
        void                    printRgf(void);
        void                    printRgf(void) const;
        double                  ssDistance(std::map<Partition*,int>& partList);

    protected:
                                Partition(void) { }
        void                    deleteSubsets(void);
        int                     flip(int x);
        Subset*                 removeElement(int x);
        void                    removeEmptySubsets(void);
        void                    setRgf(void);
        int                     numElements;
        std::set<Subset*>       subsets;
        std::vector<int>        rgf;
};

#endif
