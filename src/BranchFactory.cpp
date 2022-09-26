#include "Branch.hpp"
#include "BranchFactory.hpp"



BranchFactory::BranchFactory(void) {

}

BranchFactory::~BranchFactory(void) {

    for (std::set<Branch*>::iterator b=allocated.begin(); b != allocated.end(); b++)
        delete (*b);
}

void BranchFactory::drainPool(void) {

    for (std::vector<Branch*>::iterator b=pool.begin(); b != pool.end(); b++)
        {
        allocated.erase( *b );
        delete (*b);
        }
}

Branch* BranchFactory::getBranch(void) {

    if ( pool.empty() == true )
        {
        /* If the branch pool is empty, we allocate a new branch and return it. We
           do not need to add it to the branch pool. */
        Branch* b = new Branch;
        allocated.insert( b );
        return b;
        }
    
    // Return a branch from the branch pool, remembering to remove it from the pool.
    Branch* b = pool.back();
    pool.pop_back();
    return b;
}

void BranchFactory::returnToPool(Branch* b) {

    b->clean();
    pool.push_back( b );
}
