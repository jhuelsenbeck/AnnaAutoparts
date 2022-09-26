#include "Subset.hpp"



Subset::Subset(int n) {

    numElements = n;
    bitRep.resize(numElements);
    for (int i=0; i<numElements; i++)
        bitRep[i] = false;
}

Subset::Subset(Subset& s) {

    numElements = s.numElements;
    bitRep.resize(numElements);
    bitRep = s.bitRep;
}

Subset& Subset::operator=(Subset& s) {

    if (this != &s)
        {
        numElements = s.numElements;
        bitRep = s.bitRep;
        }
    return *this;
}

void Subset::addElement(int x) {

    bitRep[x] = true;
}

int Subset::numAssigned(void) {

    int n = 0;
    for (int i=0; i<numElements; i++)
        {
        if (bitRep[i] == true)
            n++;
        }
    return n;
}

void Subset::print(void) {

    std::cout << this << " -- ";
    int n = 0;
    for (int i=0; i<numElements; i++)
        {
        if (bitRep[i] == true)
            {
            std::cout << "1";
            n++;
            }
        else
            std::cout << "0";
        }
    std::cout << " (" << n << ")" << std::endl;
}

void Subset::removeElement(int x) {

    bitRep[x] = false;
}
