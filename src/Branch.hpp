#ifndef Branch_H
#define Branch_H

#define MIN_PROPORTION 0.000001
#define MAX_PROPORTION 0.999999

class Node;



class Branch {

    public:
                    Branch(void);
                    Branch(Node* e1, Node* e2);
                    Branch(Node* e1, Node* e2, double p);
        Branch&     operator=(const Branch& a);
        bool        operator==(const Branch& a) const;
        bool        operator<(const Branch& a) const;
        void        clean(void);
        Node*       getAncestralNode(void);
        Node*       getDescendantNode(void);
        Node*       getEnd1(void) const { return end1; }
        Node*       getEnd2(void) const { return end2; }
        double      getProportion(void) { return proportion; }
        bool        isTip(void);
        void        print(void);
        void        setEnds(Node* e1, Node* e2);
        void        setProportion(double x);

    private:
        double      proportion;
        Node*       end1;
        Node*       end2;
};

#endif
