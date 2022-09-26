#ifndef Tree_H
#define Tree_H

#include <string>
#include <sstream>
#include "Parm.h"



class Node {

	public:
                	              Node(void);  
						   Node   *getLft(void) { return lft; }
						   Node   *getRht(void) { return rht; }
						   Node   *getAnc(void) { return anc; }
						 double   getP(void) { return proportion; }
					std::string   getName(void) { return name; }
						    int   getIndex(void) { return index; }
						   bool   getIsLeaf(void) { return isLeaf; }
						    int   getFlag(void) { return flag; }
						   bool   getMarked(void) { return marked; }
						   void   setLft(Node *p) { lft = p; }
						   void   setRht(Node *p) { rht = p; }
						   void   setAnc(Node *p) { anc = p; }
						   void   setP(double x) { proportion = x; }
						   void   setName(std::string s) { name = s; }
						   void   setIndex (int x) { index = x; }
						   void   setIsLeaf(bool tf) { isLeaf = tf; }
						   void   setFlag(int x) { flag = x; }
						   void   setMarked(bool tf) { marked = tf; }
                	              
	private:
	                       Node   *lft;
						   Node   *rht;
						   Node   *anc;
					     double   proportion;
						    int   index;
					std::string   name;
						   bool   isLeaf;
						    int   flag;
						   bool   marked;

};

class Alignment;
class MbRandom;
class Tree : public Parm {

	public:
                	              Tree(MbRandom *rp, Model *mp, std::string nm, Alignment *ap, double lm, double tn);  
								  Tree(MbRandom *rp, Model *mp, std::string nm, Alignment *ap, double lm, double tn, std::string ts);
								  Tree(Tree &t);
                	              ~Tree(void);
						 double   update(void);
						   void   clone(Tree &t);
						 double   lnPriorProb(void);
					std::string   getParmString(int n);
							int   getNumNodes(void) { return numNodes; }
							int   getNumTaxa(void) { return numTaxa; }
						  Node*   getDownPassNode(int i) { return downPassSequence[i]; }
						  Node*   getRoot(void) { return root; }
					std::string   getNewick(void);
                           void   writeTree(Node *p, std::stringstream &ss);
				           void   print(void);
					std::string   getParmHeader(int n);
                	              
	private:
	                        int   numTaxa;
							int   numNodes;
					  Alignment   *alignmentPtr;
					       Node   *nodes;
						   Node   **downPassSequence;
						   Node   *root;
						   void   buildRandomTree(void);
                           void   buildTreeFromNewickDescription(std::string ts);
						    int   dex(Node *p);
						   void   getDownPassSequence(void);
						    int   getDownPassSequence(Node **p, Node *r);
						   void   passDown(Node *p, int *x);
						   void   passDown(Node *p, Node **dp, int *x);
						   void   markBranchesDown(Node *p);
						   void   showNodes(Node *p, int indent);
					     double   alpha0; // the tuning parameter for the LOCAL proposal mechanism
						 double   lambda; // the exponential parameter for the branch length
						   bool   isTreeFixed;
						 double   updateLocal(void);
						 double   updateBrlen(void);
						 double   updateTbr(void);
                            
};

#endif