#include <iostream>
#include "Branch.hpp"
#include "Node.hpp"



Branch::Branch(void) {

    end1 = NULL;
    end2 = NULL;
    proportion = 0.0;
}

Branch::Branch(Node* e1, Node* e2) {

    setEnds(e1, e2);
    proportion = 0.0;
}

Branch::Branch(Node* e1, Node* e2, double p) {

    setEnds(e1, e2);
    proportion = p;
}

Branch& Branch::operator=(const Branch& a) {

    if (this != &a)
        {
        end1 = a.end1;
        end2 = a.end2;
        proportion = a.proportion;
        }
    return *this;
}

bool Branch::operator==(const Branch& a) const {

    if ( end1 == a.end1 && end2 == a.end2 )
        return true;
    return false;
}

bool Branch::operator<(const Branch& a) const {

    if ( a.end1 < end1 )
        return true;
    else if ( a.end1 == end1 )
        {
        if ( a.end2 < end2 )
            return true;
        }
    return false;
}

void Branch::clean(void) {

    end1 = NULL;
    end2 = NULL;
    proportion = 0.0;
}

Node* Branch::getAncestralNode(void) {

    if (end1->getAncestor() == end2)
        return end2;
    else if (end2->getAncestor() == end1)
        return end1;
    return NULL;
}

Node* Branch::getDescendantNode(void) {

    if (end1->isDescendant(end2) == true)
        return end2;
    else if (end2->isDescendant(end1) == true)
        return end1;
    return NULL;
}

bool Branch::isTip(void) {

    if (end1->getIsLeaf() == true || end2->getIsLeaf() == true)
        return true;
    return false;
}

void Branch::print(void) {

    std::cout << "Branch: " << this << std::endl;
    std::cout << "   * proportion = " << proportion << std::endl;
    std::cout << "   * end1       = " << end1->getIndex() << std::endl;
    std::cout << "   * end2       = " << end2->getIndex() << std::endl;
}

void Branch::setEnds(Node* e1, Node* e2) {

    if (e1 < e2)
        {
        end1 = e1;
        end2 = e2;
        }
    else if (e2 < e1)
        {
        end1 = e2;
        end2 = e1;
        }
    else
        {
        std::cout << "Error: A branch is defined by two different nodes (" << e1 << ", " << e2 << ")" << std::endl;
        exit(1);
        }
}

void Branch::setProportion(double x) {

    if (x < MIN_PROPORTION)
        proportion = MIN_PROPORTION;
    else if (x > MAX_PROPORTION)
        proportion = MAX_PROPORTION;
    else
        proportion = x;
}
