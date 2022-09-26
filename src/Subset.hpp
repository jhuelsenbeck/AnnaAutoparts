#ifndef Subset_H
#define Subset_H

#include <iostream>
#include <set>
#include <vector>


class Subset {

    public:
                            Subset(int n);
                            Subset(Subset& s);
        Subset&             operator=(Subset& s);
        void                addElement(int x);
        std::vector<bool>&  getBitRep(void) { return bitRep; }
        bool                isElementPartOfSet(int x) { return bitRep[x]; }
        int                 numAssigned(void);
        void                print(void);
        void                removeElement(int x);
    
    protected:
                            Subset(void) { }
        int                 numElements;
        std::vector<bool>   bitRep;
};

#endif
