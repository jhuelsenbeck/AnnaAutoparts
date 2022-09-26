#include "Table.hpp"
#include "TableFactory.hpp"



TableFactory::TableFactory(void) {

}

TableFactory::~TableFactory(void) {

    for (std::set<Table*>::iterator t=allocated.begin(); t != allocated.end(); t++)
        delete (*t);
}

void TableFactory::drainPool(void) {

    for (std::vector<Table*>::iterator t=pool.begin(); t != pool.end(); t++)
        {
        allocated.erase( *t );
        delete (*t);
        }
}

Table* TableFactory::getTable(void) {

    if ( pool.empty() == true )
        {
        /* If the branch pool is empty, we allocate a new branch and return it. We
           do not need to add it to the branch pool. */
        Table* t = new Table;
        allocated.insert( t );
        return t;
        }
    
    // Return a branch from the branch pool, remembering to remove it from the pool.
    Table* t = pool.back();
    pool.pop_back();
    return t;
}

void TableFactory::returnToPool(Table* t) {

    t->clean();
    pool.push_back( t );
}
