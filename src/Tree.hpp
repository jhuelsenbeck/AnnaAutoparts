#ifndef Tree_hpp
#define Tree_hpp

#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include "NodePair.hpp"
class Branch;
class Model;
class Node;
class RandomVariable;

struct CutInfo {

    std::vector<Node*>  orphanedNodes;
    double              missingBranchLength;
    Node*               root1;
    Node*               root2;
    std::vector<Node*>  subtree1;
    std::vector<Node*>  subtree2;
    NodePair            reconnectionPair1;
    NodePair            reconnectionPair2;
};




class Tree {

    public:
                                                    Tree(void) = delete;
                                                    Tree(Tree& t);
                                                    Tree(Model* m, std::vector<std::string> tn);
                                                    Tree(Model* m, std::vector<std::string> tn, std::string newickStr);
                                                   ~Tree(void);
        Tree&                                       operator=(Tree& t);
        Branch*                                     addBranch(Node* e1, Node* e2);
        Branch*                                     addBranch(Node* e1, Node* e2, double x);
        Node*                                       addNode(void);
        double                                      branchProportionSum(void);
        void                                        chooseBackbone(RandomVariable* rng, std::vector<Node*>& backboneNodes);
        NodePair                                    choosePair(std::map<NodePair,double,CompNodePair>& vals);
        void                                        calculateReconnectionProbabilities(CutInfo& info, std::map<NodePair,double,CompNodePair>& vals, double heat);
        Branch*                                     findBranch(Node* e1, Node* e2);
        Node*                                       findNodeIndexed(int idx);
        NodePair                                    findOriginalPair(std::map<NodePair,double,CompNodePair>& vals, CutInfo& info);
        void                                        flipAllActiveCls(void);
        void                                        flipAllActiveCls(Node* p);
        void                                        flipAllActiveTps(void);
        std::map<NodePair,Branch*,CompNodePair>&    getBranches(void) { return branches; }
        std::vector<Node*>&                         getDownPassSequence(void) { return downPassSequence; }
        int                                         getNumLeaves(void) { return numLeaves; }
        std::string                                 getNewick(void);
        int                                         getNumNodes(void) { return (int)nodes.size(); }
        int                                         getNumTaxa(void) { return (int)taxonNames.size(); }
        Node*                                       getRoot(void) { return root; }
        std::vector<std::string>                    getTaxonNames(void) { return taxonNames; }
        void                                        initalizeDownPassSequence(void);
        bool                                        isBinary(void);
        void                                        normalizeBranchProportions(void);
        void                                        print(void);
        void                                        print(std::string h);
        void                                        printBranches(std::string header);
        Branch*                                     randomBranch(RandomVariable* rng);
        Branch*                                     randomInteriorBranch(RandomVariable* rng);
        double                                      randomlyCutTree(RandomVariable* rng, CutInfo& info);
        double                                      reconnect(Node* p1, Node* p2, CutInfo& info);
        void                                        removeAllBranches(void);
        void                                        removeBranch(Node* e1, Node* e2);
        void                                        rootOnTaxon(int idx);
        void                                        setAllClUpdateFlags(bool tf);
        void                                        setAllClUpdateFlags(bool tf, Node* p);
        void                                        setAllTpUpdateFlags(bool tf);
        void                                        setAllTpUpdateFlags(bool tf, Node* p);
        void                                        setNumLeaves(int x) { numLeaves = x; }
        void                                        setRoot(Node* p) { root = p; }
        void                                        setTaxonNames(std::vector<std::string> tn) { taxonNames = tn; }
        void                                        setModel(Model* m) {modelPtr = m;}
    
    private:
        void                                        clone(Tree& t);
        void                                        deleteAllNodes(void);
        std::vector<std::string>                    parseNewickString(std::string ns);
        void                                        passDown(Node* p, Node* from);
        void                                        passDown(Node* p, Node* from, std::vector<Node*>& dp);
        void                                        printNode(Node* p, int indent);
        void                                        tipDescendantsFrom(Node* p, Node* from, std::set<int>& tips);
        void                                        writeTree(Node* p, std::stringstream& ss);
        std::vector<Node*>                          nodes;
        std::map<NodePair,Branch*,CompNodePair>     branches;
        NodePair                                    branchKey;
        std::vector<std::string>                    taxonNames;
        std::vector<Node*>                          downPassSequence;
        Node*                                       root;
        int                                         numLeaves;
        Model*                                      modelPtr;
};

#endif
