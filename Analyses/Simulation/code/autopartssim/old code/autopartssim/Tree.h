#ifndef Tree_H
#define Tree_H

#include <string>
#include <vector>


class MbRandom;
class Node;
class Tree {

	public:
                                        Tree(int nt, MbRandom* rp);
                                       ~Tree(void);
                                 void   getDownPassSequence(void);
                          std::string   getNewick(void);
                                Node*   getRoot(void) { return root; }
                  std::vector<Node*>&   getTraversalSequence(void) { return downPassSequence; }
                                 void   print(void);

	private:
                                 void   buildRandomTree(void);
                                  int   dex(Node* p);
                                 void   passDn(Node* p);
                                 void   writeTree(Node* p, std::stringstream &ss);
                                Node*   nodes;
                   std::vector<Node*>   downPassSequence;
                                Node*   root;
                            MbRandom*   ranPtr;
                                  int   numTaxa;
                                  int   numNodes;
};

#endif