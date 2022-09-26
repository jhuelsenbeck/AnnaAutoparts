#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include "MbRandom.h"
#include "Node.h"
#include "Tree.h"



Tree::Tree(int nt, MbRandom* rp) {

    numTaxa  = nt;
    ranPtr   = rp;
    numNodes = 2 * numTaxa - 2;
    
    buildRandomTree();
}

Tree::~Tree(void) {

    delete [] nodes;
}

void Tree::buildRandomTree(void) {

    // allocate the nodes
    nodes = new Node[numNodes];

    // set the index uniquely for each node
    for (int i=0; i<numNodes; i++)
        nodes[i].setIndex(i);
        
    // build a randomly bifurcating unrooted tree
    Node* p = &nodes[0];
    root = p;
    Node* q = &nodes[numTaxa];
    p->setLft(q);
    q->setAnc(p);
    p = q;
    q = &nodes[1];
    p->setLft(q);
    q->setAnc(p);
    q = &nodes[2];
    p->setRht(q);
    q->setAnc(p);
    std::vector<Node*> availableNodes;
    availableNodes.push_back( &nodes[1] );
    availableNodes.push_back( &nodes[2] );
    availableNodes.push_back( &nodes[numTaxa] );
    int nextInteriorNode = numTaxa+1;
        
    for (int i=3; i<numTaxa; i++)
        {
        // pick a branch at random
        int whichNode = (int)( ranPtr->uniformRv()*availableNodes.size() );
        p = availableNodes[whichNode];
        Node* a = p->getAnc();
        q = &nodes[nextInteriorNode];
        nextInteriorNode++;
        Node* r = &nodes[i];
        if ( p->getAnc()->getLft() == p )
            {
            // p leans to the left
            a->setLft(q);
            q->setAnc(a);
            q->setLft(p);
            p->setAnc(q);
            q->setRht(r);
            r->setAnc(q);
            }
        else 
            {
            // p leans to the right
            a->setRht(q);
            q->setAnc(a);
            q->setRht(p);
            p->setAnc(q);
            q->setLft(r);
            r->setAnc(q);
            }
        availableNodes.push_back( q );
        availableNodes.push_back( r );
        }
 
    // get the traversal sequence
    getDownPassSequence();
    
    // initialize the branch lengths by drawing from an exponential repeatedly
    double sum = 0.0;
    for (int i=0; i<downPassSequence.size(); i++)
        {
        Node* p = downPassSequence[i];
        if (p != root)
            {
            p->setProportion( ranPtr->exponentialRv(1.0) );
            sum += p->getProportion();
            }
        }
    for (int i=0; i<downPassSequence.size(); i++)
        {
        Node* p = downPassSequence[i];
        if (p != root)
            p->setProportion( p->getProportion() / sum );
        }
}

int Tree::dex(Node* p) {

    if (p == NULL)
        return -1;
    return p->getIndex();
}

void Tree::getDownPassSequence(void) {

    downPassSequence.clear();
    passDn(root);
}

std::string Tree::getNewick(void) {

	std::stringstream ss;
	writeTree(root->getLft(), ss);
	std::string newick = ss.str();
	return newick;
}

void Tree::print(void) {

    for (int i=0; i<numNodes; i++)
        {
        Node* p = &nodes[i];
        std::cout << i << " -- " << p->getIndex() << " (" << dex(p->getLft()) << " " << dex(p->getRht()) << " " << dex(p->getAnc()) << ")";
        std::cout << " " << std::fixed << std::setprecision(4) << p->getProportion();
        if (p == root)
            std::cout << " <- Root" << std::endl;
        else 
            std::cout << std::endl;
        }
}

void Tree::writeTree(Node* p, std::stringstream &ss) {

	if (p != NULL)
		{
		
		if (p->getLft() == NULL && p->getRht() == NULL)
			{
			ss << p->getIndex() << ":" << std::fixed << std::setprecision(6) << p->getProportion();
			}
		else
			{
			if (p->getAnc() != NULL)
				{
				ss << "(";
				}
			writeTree(p->getLft(), ss);
			ss << ",";
			writeTree(p->getRht(), ss);	
			if (p->getAnc() != NULL)
				{
				if (p->getAnc()->getAnc() == NULL)
					{
					ss << "," << p->getAnc()->getIndex() << ":" << std::fixed << std::setprecision(6) << p->getProportion();
					}
				
				if (p->getAnc()->getAnc() != NULL)
					ss << "):" << std::fixed << std::setprecision(6) << p->getProportion();
				else
					ss << ")";
				}
			}
		}

}

void Tree::passDn(Node* p) {

    if (p != NULL)
        {
        passDn(p->getLft());
        passDn(p->getRht());
        downPassSequence.push_back( p );
        }
}
