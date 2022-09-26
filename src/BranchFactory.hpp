#ifndef BranchFactory_H
#define BranchFactory_H

#include <set>
#include <vector>
class Branch;



class BranchFactory {

    public:
        static BranchFactory&   branchFactoryInstance(void)
                                    {
                                    static BranchFactory singleBranchFactory;
                                    return singleBranchFactory;
                                    }
        void                    drainPool(void);
        Branch*                 getBranch(void);
        int                     getNumAllocated(void) { return (int) allocated.size(); }
        int                     getNumInPool(void) { return (int)pool.size(); }
        void                    returnToPool(Branch* b);

    private:
                                BranchFactory(void);
                                BranchFactory(const BranchFactory&);
                                BranchFactory& operator=(const BranchFactory&);
                               ~BranchFactory(void);
        std::vector<Branch*>    pool;
        std::set<Branch*>       allocated;
};

#endif
