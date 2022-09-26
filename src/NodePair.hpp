#ifndef NodePair_hpp
#define NodePair_hpp

class Node;



class NodePair {

    public:
                NodePair(void);
                NodePair(Node* n1, Node* n2);
                NodePair(const Node* n1, const Node* n2);
        Node*   getEnd1(void) { return e1; }
        Node*   getEnd2(void) { return e2; }
        Node*   getEnd1(void) const { return e1; }
        Node*   getEnd2(void) const { return e2; }
        void    print(void);
        void    setEnds(Node* n1, Node* n2);
        
    private:
        Node*   e1;
        Node*   e2;
};

struct CompNodePair {

    bool operator()(const NodePair& np1, const NodePair& np2) const {
        
        if ( np1.getEnd1() < np2.getEnd1() )
            return true;
        else if ( np1.getEnd1() == np2.getEnd1() )
            {
            if ( np1.getEnd2() < np2.getEnd2() )
                return true;
            }
        return false;
        }
        
    bool operator()(NodePair& np1, NodePair& np2) const {
        
        if ( np1.getEnd1() < np2.getEnd1() )
            return true;
        else if ( np1.getEnd1() == np2.getEnd1() )
            {
            if ( np1.getEnd2() < np2.getEnd2() )
                return true;
            }
        return false;
        }
};

#endif
