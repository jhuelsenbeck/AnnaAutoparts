#ifndef TableFactory_hpp
#define TableFactory_hpp

#include <set>
#include <vector>
class Table;



class TableFactory {

    public:
        static TableFactory&    tableFactory(void)
                                    {
                                    static TableFactory singleTableFactory;
                                    return singleTableFactory;
                                    }
        void                    drainPool(void);
        Table *                 getTable(void);
        int                     getNumAllocated(void) { return (int) allocated.size(); }
        int                     getNumInPool(void) { return (int)pool.size(); }
        void                    returnToPool(Table* t);

    private:
                                TableFactory(void);
                                TableFactory(const TableFactory& tf) = delete;
                                TableFactory& operator=(const TableFactory&) = delete;
                               ~TableFactory(void);
        std::vector<Table*>     pool;
        std::set<Table*>        allocated;
};

#endif
