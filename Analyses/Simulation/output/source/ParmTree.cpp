#include <iostream>
#include <iomanip>
#include <string>
#include "Alignment.h"
#include "MbRandom.h"
#include "Model.h"
#include "ParmTree.h"
#include "StateSets.h"

#define	BRLENS_MIN					0.000001
#define	BRLENS_MAX					10.0
#define	BRPROP_MIN					0.000001
#define	BRPROP_MAX					10.0



Node::Node(void) {

	lft        = NULL;
	rht        = NULL;
	anc        = NULL;
	proportion = 0.0;
	index      = 0;
	name       = "";
	isLeaf     = false;
	flag       = 0;
	marked     = false;
}

Tree::Tree(MbRandom *rp, Model *mp, std::string nm, Alignment *ap, double lm, double tn) : Parm(rp, mp, nm) {

	isTreeFixed  = false;
	alpha0       = tn;
	lambda       = lm;
	alignmentPtr = ap;
	numTaxa      = alignmentPtr->getNumTaxa();
	numNodes     = 2 * numTaxa - 2;
	buildRandomTree();
}

Tree::Tree(MbRandom *rp, Model *mp, std::string nm, Alignment *ap, double lm, double tn, std::string ts) : Parm(rp, mp, nm) {

	isTreeFixed  = true;
	alpha0       = tn;
	lambda       = lm;
	alignmentPtr = ap;
	numTaxa      = alignmentPtr->getNumTaxa();
	numNodes     = 2 * numTaxa - 2;
	buildTreeFromNewickDescription(ts);
}

Tree::Tree(Tree &t) : Parm(t.ranPtr, t.modelPtr, t.parmName) {

	isTreeFixed      = false;
	alpha0           = 0.0;
	lambda           = 0.0;
	alignmentPtr     = NULL;
	numTaxa          = 0;
	numNodes         = 0;
	nodes            = NULL;
	downPassSequence = NULL;
	root             = NULL;
	clone(t);
}

Tree::~Tree(void) {

	delete [] nodes;
	delete [] downPassSequence;
}

void Tree::buildRandomTree(void) {

	/* initialize */
	int nextTip = 0;
	int nextInt = numTaxa;
	nodes = new Node[numNodes];
	downPassSequence = new Node*[numNodes];
	int numAvailableNodes = 0;
	Node **availableNodes = new Node*[numNodes];
	for (int i=0; i<numNodes; i++)
		nodes[i].setIndex( i );
	for (int i=0; i<numTaxa; i++)
		{
		nodes[i].setIsLeaf( true );
		nodes[i].setName( alignmentPtr->getTaxonName(i) );
		}
	
	/* build three-species tree */
	Node *p = &nodes[nextTip++];
	Node *q = &nodes[nextInt++];
	availableNodes[numAvailableNodes++] = q;
	root = p;
	p->setLft( q );
	q->setAnc( p );
	p = q;
	q = &nodes[nextTip++];
	availableNodes[numAvailableNodes++] = q;
	p->setLft( q );
	q->setAnc( p );
	q = &nodes[nextTip++];
	availableNodes[numAvailableNodes++] = q;
	p->setRht( q );
	q->setAnc( p );
	
	/* build the remaining portion of  the tree */
	for (int i=4; i<=numTaxa; i++)
		{
		/* pick a node at random */
		int whichNode = (int)(ranPtr->uniformRv()*numAvailableNodes);
		p = availableNodes[whichNode];
		Node *pAnc = p->getAnc();
		
		Node *newTip = &nodes[nextTip++];
		Node *newInt = &nodes[nextInt++];
		
		if (pAnc->getLft() == p)
			{
			pAnc->setLft( newInt );
			newInt->setAnc( pAnc );
			newInt->setLft( p );
			p->setAnc( newInt );
			newInt->setRht( newTip );
			newTip->setAnc( newInt );
			}
		else
			{
			pAnc->setRht( newInt );
			newInt->setAnc( pAnc );
			newInt->setRht( p );
			p->setAnc( newInt );
			newInt->setLft( newTip );
			newTip->setAnc( newInt );
			}
			
		availableNodes[numAvailableNodes++] = newInt;
		availableNodes[numAvailableNodes++] = newTip;
			
		}
		
	/* initialize the branch-length proportions */
	double sum = 0.0;
	for (int i=0; i<numNodes; i++)
		{
		p = &nodes[i];
		if (p->getAnc() != NULL)
			{
			double x = ranPtr->exponentialRv(lambda);
			sum += x;
			p->setP( x );
			}
		}
	for (int i=0; i<numNodes; i++)
		{
		p = &nodes[i];
		if (p->getAnc() != NULL)
			p->setP( p->getP()/sum );
		}
		
	/* remember post-order traversal sequence */
	getDownPassSequence();
	
	/* free memory */
	delete [] availableNodes;
	
#	if 0
	print();
#	endif

}

void Tree::buildTreeFromNewickDescription(std::string ts) {

	/* parse the tree string, and put each token into a vector of strings */
	std::vector<std::string> parsedNewick;
	std::string temp = "";
	bool readingBrlen = false;
	int nt = 0;
	for (int i=0; i<ts.size(); i++)
		{
		char c = ts[i];
		if ( c == ' ' )
			continue;
		if ( c == '(' || c == ')' || c == ',' || c == ':' || c == ';' )
			{
			temp = c;
			parsedNewick.push_back( temp );
			if ( c == ':' )
				readingBrlen = true;
			else
				readingBrlen = false;
			}
		else
			{
			/* the character is part of a taxon name */
			int j = i;
			std::string taxonName = "";
			while ( ts[j] != '(' && ts[j] != ')' && ts[j] != ',' && ts[j] != ':' && ts[j] != ';' )
				{
				taxonName += ts[j];
				j++;
				}
			parsedNewick.push_back( taxonName );
			i = j - 1;
			if ( readingBrlen == false )
				nt++;
			readingBrlen = false;
			}
		if ( c == ';' )
			break;
		}

	/* check that the number of taxa in the tree description is the same as the
	   number of taxa in the alignment */
	if ( nt != numTaxa )
		{
		std::cerr << "ERROR: The tree file is not the right size" << std::endl;
		std::cout << "nt = " << nt << " numTaxa = " << numTaxa << std::endl;
		exit(1);
		}
		
	/* allocate the nodes and down pass sequence */
	nodes = new Node[numNodes];
	if ( !nodes )
		{
		std::cerr << "Problem allocating nodes!" << std::endl;
		exit(1);
		}
	downPassSequence = new Node*[numNodes];
	if ( !downPassSequence )
		{
		std::cerr << "Problem allocating downPass!" << std::endl;
		exit(1);
		}

	/* set the indices for all of the nodes */
	for (int i=0; i<numNodes; i++)
		nodes[i].setIndex( i );
	
	/* build up the tree, using the information stored in the parsed vector */
	int nextInteriorNode = numTaxa;
	Node *p;
	int n = 0;
	double inputTreeLength = 0.0;
	for (std::vector<std::string>::iterator t=parsedNewick.begin(); t != parsedNewick.end(); t++)
		{
		//cout << (*t) << endl;
		if ( (*t) == "(" )
			{
			/* add a new interior node */
			if (n == 0)
				{
				p = &nodes[nextInteriorNode++];
				}
			else
				{
				Node *q = &nodes[nextInteriorNode++];
				if (p->getLft() == NULL)
					{
					p->setLft( q );
					q->setAnc( p );
					p = q;
					}
				else if (p->getRht() == NULL)
					{
					p->setRht( q );
					q->setAnc( p );
					p = q;
					}
				else if (p->getAnc() == NULL)
					{
					p->setAnc( q );
					q->setLft( p );
					p = q;
					p->setFlag( true );
					}
				else
					{
					std::cout << "ERROR: Problem reading the Newick-formatted tree (1)" << std::endl;
					exit(1);
					}
				}
			n++;
			readingBrlen = false;
			}
		else if ( (*t) == ")" )
			{
			/* we hit a right parantheses, so we should go down the tree...unless
			   we should go up! */
			if (p->getFlag() == false && p->getAnc() != NULL)
				p = p->getAnc();
			else if (p->getFlag() == true && p->getLft() != NULL)
				p = p->getLft();
			else
				{
				std::cout << "ERROR: Problem reading the Newick-formatted tree (2)" << std::endl;
				exit(1);
				}
			readingBrlen = false;
			}
		else if ( (*t) == "," )
			{
			/* we hit a comma, so we should go down the tree...unless
			   we should go up! */
			if (p->getFlag() == false && p->getAnc() != NULL)
				p = p->getAnc();
			else if (p->getFlag() == true && p->getLft() != NULL)
				p = p->getLft();
			else
				{
				std::cout << "ERROR: Problem reading the Newick-formatted tree (3)" << std::endl;
				exit(1);
				}
			readingBrlen = false;
			}
		else if ( (*t) == ":" )
			{
			readingBrlen = true;
			}
		else if ( (*t) == ";" )
			{
			/* We are at the end of the tree description. I guess we don't have
			   to do anything. */
			}
		else
			{
			if (readingBrlen == false)
				{
				/* read in a taxon name, and add the node to the tree */
				int theIndex;
				std::istringstream buf(*t);
				buf >> theIndex;
				theIndex--;
				std::string theName = alignmentPtr->getTaxonName(theIndex);
				
				//string theName = (*t);
				//int theIndex = alignmentPtr->getTaxonIndex( theName );
				
				Node *q = &nodes[theIndex];
				q->setName( theName );
				q->setIsLeaf( true );
				
				if (p->getLft() == NULL)
					{
					p->setLft( q );
					q->setAnc( p );
					p = q;
					}
				else if (p->getRht() == NULL)
					{
					p->setRht( q );
					q->setAnc( p );
					p = q;
					}
				else if (p->getAnc() == NULL)
					{
					p->setAnc( q );
					q->setLft( p );
					p = q;
					p->setFlag( true );
					}
				else
					{
					std::cout << "ERROR: Problem reading the Newick-formatted tree (4)" << std::endl;
					exit(1);
					}
				n++;
				}
			else
				{
				/* reading a branch length */
				double x = 0.0;
				std::istringstream buf(*t);
				buf >> x;
				if (x < 0.00001)
					x = 0.0001;
				inputTreeLength +=  x;
				
				if (p->getFlag() == false)
					p->setP( x );
				else if (p->getFlag() == true && p->getLft() != NULL)
					p->getLft()->setP( x );
				else
					{
					std::cout << "ERROR: I have no clue where to put this branch length" << std::endl;
					exit(1);
					}

				
				readingBrlen = false;
				}
			}
		}

	/* set the pointer to the root */
	Node *q = p;
	while (q->getAnc() != NULL)
		q = q->getAnc();
	root = q;
	
	/* get the post-order traversal sequence */
	getDownPassSequence();
	
	/* initialize the branch length proportions */
	double sum = 0.0;
	for (int i=0; i<numNodes; i++)
		{
		p = getDownPassNode( i );
		if (p->getAnc() != NULL)
			sum += p->getP();
		}
	for (int i=0; i<numNodes; i++)
		{
		p = &nodes[i];
		if (p->getAnc() != NULL)
			p->setP( p->getP()/sum );
		}

#	if 0
	print();
	cout << "T = " << inputTreeLength << endl;
#	endif
}

double Tree::lnPriorProb(void) {

	std::vector<double> brAlpha(2*getNumTaxa() - 3);
	std::vector<double> brProps(2*getNumTaxa() - 3);
	for (int n=0, k=0; n<numNodes; n++)
		{
		Node *p = getDownPassNode(n);
		if (p->getAnc() != NULL)
			{
			brAlpha[k  ] = 1.0;
			brProps[k++] = p->getP();
			}
		}
	return ranPtr->lnDirichletPdf( brAlpha, brProps );
}

std::string Tree::getParmString(int n) {

	std::string tempString = "";
	char tempCharStr[10];
	sprintf(tempCharStr, "%d", n);
	std::string anotherTempStr = tempCharStr;
	return tempString = "tree sample_" + anotherTempStr + " = " + getNewick() + ";";
}

std::string Tree::getParmHeader(int n) {

	std::string tempString = "";
	return tempString;
}

std::string Tree::getNewick(void) {

	std::stringstream ss;
	writeTree(root->getLft(), ss);
	std::string newick = ss.str();
	return newick;
}

void Tree::writeTree(Node *p, std::stringstream &ss) {

	if (p != NULL)
		{
		
		if (p->getLft() == NULL && p->getRht() == NULL)
			{
			ss << p->getName() << ":" << std::fixed << std::setprecision(6) << p->getP();
			}
		else
			{
			if (p->getAnc() != NULL)
				{
				ss << "(";
				}
			writeTree (p->getLft(), ss);
			ss << ",";
			writeTree (p->getRht(), ss);	
			if (p->getAnc() != NULL)
				{
				if (p->getAnc()->getAnc() == NULL)
					{
					ss << "," << p->getAnc()->getName() << ":" << std::fixed << std::setprecision(6) << p->getP();
					}
				
				if (p->getAnc()->getAnc() != NULL)
					ss << "):" << std::fixed << std::setprecision(6) << p->getP();
				else
					ss << ")";
				}
			}
		}

}

double Tree::update(void) {

	double lnP = 0.0;
	if (isTreeFixed == false)
		{
		if (ranPtr->uniformRv() < 0.2)
			lnP = updateLocal();
		else
			lnP = updateTbr();
		}
	else
		{
		lnP = updateBrlen();
		}
	return lnP;
}

double Tree::updateBrlen(void) {

	/* randomly pick a branch */
	Node *u;
	do
		{
		u = downPassSequence[(int)(ranPtr->uniformRv()*numNodes)];
		} while ( u->getAnc() == NULL );
	
	/* select new proportions */
	double oldP = u->getP();
	std::vector<double> alp(2);
	std::vector<double> z(2);
	alp[0] = oldP * alpha0;
	alp[1] = (1.0 - oldP) * alpha0;
	//cout << fixed << setprecision(25) << alp[0] << " " << alp[1] << endl;
	ranPtr->dirichletRv(alp, z);
	double newP = z[0];
	if (newP < 0.000001)
		{
		newP = 0.000001;
		z[0] = newP;
		z[1] = 1.0 - newP;
		}
	//cout << fixed << setprecision(25) << newP << " " << 1.0 - newP << endl;
	double lnForwardProb = ranPtr->lnDirichletPdf(alp, z);
	alp[0] = newP * alpha0;
	alp[1] = (1.0 - newP) * alpha0;
	z[0] = oldP;
	z[1] = 1.0 - oldP;
	double lnReverseProb = ranPtr->lnDirichletPdf(alp, z);
	
	/* update the branch lengths */
	double sum = 0.0;
	for (int n=0; n<numNodes; n++)
		{
		Node *p = downPassSequence[n];
		if (p->getAnc() != NULL)
			{
			double v = p->getP();
			if ( p == u )
				p->setP( newP );
			else
				p->setP( v*((1.0-newP)/(1.0-oldP)) );
			sum += p->getP();
			}
		}
	for (int n=0; n<numNodes; n++) // rescale branch lengths, just in case (not much should be going on here)
		{
		Node *p = downPassSequence[n];
		if (p != root)
			p->setP( p->getP()/sum );
		}
	
	return (lnReverseProb - lnForwardProb) + (numNodes-2)*log((1.0-newP)/(1.0-oldP));
}

double Tree::updateLocal(void) {

#	if 1

	/* pick an internal branch at random */
	Node *u;
	do
		{
		u = downPassSequence[(int)(ranPtr->uniformRv()*numNodes)];
		} while( !(u->getIsLeaf() == false && root->getLft() != u) );
	Node *v = u->getAnc();
	bool isLeft = false;
	if (v->getLft() == u)
		isLeft = true;
		
	/* pick two other branches incident to either end of branch a, to form the backbone of the move */
	Node *a, *b, *c, *d;
	if (ranPtr->uniformRv() < 0.5)
		{
		a = u->getLft();
		c = u->getRht();
		}
	else
		{
		a = u->getRht();
		c = u->getLft();
		}
	bool isBelow = false;
	if (ranPtr->uniformRv() < 0.5)
		{
		b = v->getAnc();
		if (isLeft == true)
			d = v->getRht();
		else
			d = v->getLft();
		isBelow = true;
		}
	else
		{
		d = v->getAnc();
		if (isLeft == true)
			b = v->getRht();
		else
			b = v->getLft();
		}
	
	/* store path in vector */
	Node *path[3];
	path[0] = a;
	path[1] = u;
	if (isBelow == true)
		path[2] = v;
	else
		path[2] = b;
		
	/* pick a new path length */
	double oldM = a->getP() + u->getP();
	if (isBelow == true)
		oldM += v->getP();
	else
		oldM += b->getP();
	std::vector<double> alp(2);
	std::vector<double> z(2);
	alp[0] = oldM * alpha0;
	alp[1] = (1.0 - oldM) * alpha0;
	ranPtr->dirichletRv(alp, z);
	double newM = z[0];
	double lnForwardProb = ranPtr->lnDirichletPdf(alp, z);
	alp[0] = newM * alpha0;
	alp[1] = (1.0 - newM) * alpha0;
	z[0] = oldM;
	z[1] = 1.0 - oldM;
	double lnReverseProb = ranPtr->lnDirichletPdf(alp, z);
	
	/* calculate the log of the Hastings ratio and Jacobian here */
	double lnProposalProb = lnReverseProb - lnForwardProb;
	lnProposalProb += 2.0 * (log(newM) - log(oldM)) + (numNodes-4) * (log(1.0-newM) - log(1.0-oldM));
	
	/* reset all of the branch lengths */
	double sum = 0.0;
	for (int n=0; n<numNodes; n++)
		{
		Node *p = downPassSequence[n];
		if (p != root)
			{
			double v = p->getP();
			if ( p == path[0] || p == path[1] || p == path[2] )
				p->setP( v*(newM/oldM) );
			else
				p->setP( v*((1.0-newM)/(1.0-oldM)) );
			sum += p->getP();
			}
		}
	for (int n=0; n<numNodes; n++) // rescale branch lengths, just in case (not much should be going on here)
		{
		Node *p = downPassSequence[n];
		if (p != root)
			p->setP( p->getP()/sum );
		}
		
	/* randomly pick one of the two branches that are incident to the path */
	Node *nodeToMove = c;
	if (ranPtr->uniformRv() < 0.5)
		nodeToMove = d;

	/* randomly reattach the selected node to the path */
	double newPos = ranPtr->uniformRv() * newM;
	bool topologyChanged = false;
	if (nodeToMove == c)
		{
		/* node c is moved */
		double criticalVal = a->getP() + u->getP();
		if (newPos < criticalVal)
			{
			/* no topology change */
			a->setP(newPos);
			u->setP(criticalVal-newPos);
			}
		else
			{
			/* topology change */
			a->setP( a->getP()+u->getP() );
			if (v->getLft() == u)
				v->setLft(a);
			else
				v->setRht(a);
			a->setAnc(v);
			if (isBelow == true)
				{
				if (b->getLft() == v)
					{
					b->setLft(u);
					u->setAnc(b);
					u->setLft(v);
					v->setAnc(u);
					u->setRht(c);
					c->setAnc(u);
					v->setP(newPos - criticalVal);
					u->setP(newM - newPos);
					}
				else
					{
					b->setRht(u);
					u->setAnc(b);
					u->setLft(c);
					c->setAnc(u);
					u->setRht(v);
					v->setAnc(u);
					}
				v->setP(newPos - criticalVal);
				u->setP(newM - newPos);
				}
			else
				{
				if (v->getLft() == b)
					{
					v->setLft(u);
					u->setAnc(v);
					u->setLft(b);
					b->setAnc(u);
					u->setRht(c);
					c->setAnc(u);
					}
				else
					{
					v->setRht(u);
					u->setAnc(v);
					u->setLft(c);
					c->setAnc(u);
					u->setRht(b);
					b->setAnc(u);
					}
				b->setP(newPos - criticalVal);
				u->setP(newM - newPos);
				}
			topologyChanged = true;
			}
		}
	else
		{
		/* node d is moved */
		double criticalVal = a->getP();
		if (newPos > criticalVal)
			{
			/* no topology change */
			if (isBelow == true)
				{
				u->setP(newPos - criticalVal);
				v->setP(newM - newPos);
				}
			else
				{
				u->setP(newPos - criticalVal);
				b->setP(newM - newPos);
				}
			}
		else
			{
			/* topology change */
			if (isBelow == true)
				{
				u->setP( u->getP()+v->getP() );
				if (b->getLft() == v)
					{
					b->setLft(u);
					u->setAnc(b);
					}
				else
					{
					b->setRht(u);
					u->setAnc(b);
					}
				if (u->getLft() == a)
					{
					u->setLft(v);
					v->setAnc(u);
					v->setLft(a);
					a->setAnc(v);
					v->setRht(d);
					d->setAnc(v);
					}
				else
					{
					u->setRht(v);
					v->setAnc(u);
					v->setLft(d);
					d->setAnc(v);
					v->setRht(a);
					a->setAnc(v);
					}
				a->setP(newPos);
				v->setP(criticalVal - newPos);
				}
			else
				{
				b->setP( b->getP()+u->getP() );
				v->setLft(u);
				v->setRht(a);
				u->setAnc(v);
				a->setAnc(v);
				u->setLft(b);
				u->setRht(c);
				b->setAnc(u);
				c->setAnc(u);
				u->setP(criticalVal - newPos);
				a->setP(newPos);
				}			
			topologyChanged = true;
			}
		}
	
	if (topologyChanged == true)
		getDownPassSequence();
	
	return lnProposalProb;

#	else
	/* we need to keep track of the forward and reverse proposal probabilities so we can
	   calculate the hastings ratio */
	double lnProposalProb;

	bool topologyHasChanged = false;
	
	/* pick an internal branch */
	Node *v;
	bool goodNode = false;
	do
		{
		v = downPassSequence[ (int)(ranPtr->uniformRv()*numNodes) ];
		if ( v->getLft() != NULL && v->getRht() != NULL && v->getAnc() != NULL )
			{
			if (v->getAnc() != root)
				goodNode = true;
			}
		} while( goodNode == false );
		
	/* set up pointers for crown part */
	Node *c, *d;
	if (ranPtr->uniformRv() < 0.5)
		{
		c = v->getLft();
		d = v->getRht();
		}
	else
		{
		c = v->getRht();
		d = v->getLft();
		}

	/* set up pointers for root part */
	Node *u = v->getAnc();
	bool directionUp;
	Node *a, *b;
	if (ranPtr->uniformRv() < 0.5)
		{
		directionUp = true;
		if (u->getLft() == v)
			a = u->getRht();
		else
			a = u->getLft();
		b = u->getAnc();
		}
	else
		{
		directionUp = false;
		if (u->getLft() == v)
			b = u->getRht();
		else
			b = u->getLft();
		a = u->getAnc();
		}

	/* find length multiplication factor */
	double lenFactor = exp(tuning * (ranPtr->uniformRv() - 0.5));

	/* store old and new path length as well as old x and y */
	double oldM = c->getLength() + v->getLength();
	double x = 0.0;
	if (directionUp == true)
		{
		oldM += a->getLength();
		x = a->getLength();
		}
	else
		{
		oldM += u->getLength();
		x = u->getLength();
		}

	double y = x + v->getLength();
	double newM = oldM * lenFactor;
	
	/* update the transition probability update flags */
	/*c->setUpdateTi( true );
	c->flipActiveTi();
	v->setUpdateTi( true );
	v->flipActiveTi();
	if (directionUp == true)
		{
		a->setUpdateTi( true );
		a->flipActiveTi();
		}
	else
		{
		u->setUpdateTi( true );
		u->flipActiveTi();
		}*/
	
	/* adjust proposal and prior ratio based on length modification */
	/* and insertion mechanism */
	lnProposalProb = 3.0 * log(lenFactor);

	/* pick dangly to move and new attachment point */
	if (ranPtr->uniformRv() < 0.5)
		{
		/* choose new x */
		x = ranPtr->uniformRv() * newM;
		y *= lenFactor;
		}
	else
		{
		/* choose new y */
		y = ranPtr->uniformRv() * newM;
		x *= lenFactor;
		}

	/* make topology move if necessary and then set branch lengths */
	if (x > y)
		{
		/* topology has changed */
		topologyHasChanged = true;
		/* detach v and d */
		/* this scheme differs from that used by Larget and Simon but is more
		   convenient because it avoids tree rotations */
		if (u->getLft() == v)
			u->setLft( c );
		else
			u->setRht( c );
		c->setAnc( u );
		if (directionUp == true)
			{
			/* place v and d below a */
			if (v->getLft() == d)
				v->setRht( a );
			else
				v->setLft( a );
			a->setAnc( v );
			if (u->getLft() == a)
				u->setLft( v );
			else
				u->setRht( v );
			/* v->anc is already u */
			/* adjust lengths */
			c->setLength( newM - x );
			v->setLength( x - y );
			a->setLength( y );
			}
		else
			{
			/* place v and d below u */
			if (v->getLft() == d)
				v->setRht( u );
			else
				v->setLft( u );
			u->setAnc( v );
			v->setAnc( a );
			if (a->getLft() == u)
				a->setLft( v );
			else
				a->setRht( v );
			/* adjust lengths */
			c->setLength( newM - x );
			u->setLength( x - y );
			v->setLength( y );
			}
		}
	else
		{
		/* topology has not changed */
		c->setLength( newM - y );
		v->setLength( y - x );
		if (directionUp == true)
			a->setLength( x );
		else
			u->setLength( x );
		}

	/* update the conditional likelihood update flags */
	/*Node *p = c->getAnc();
	while ( p != root )
		{
		p->setUpdateCl( true );
		p->flipActiveCl();
		p = p->getAnc();
		}*/
		
	/* check branch lengths */
	double minV = BRLENS_MIN;
	double maxV = BRLENS_MAX;
	if (c->getLength() > maxV)
		c->setLength( maxV );
	if (v->getLength() > maxV)
		v->setLength( maxV );
	if (c->getLength() < minV)
		c->setLength( minV );
	if (v->getLength() < minV)
		v->setLength( minV );
	if (directionUp == true)
		{
		if (a->getLength() > maxV)
			a->setLength( maxV );
		if (a->getLength() < minV)
			a->setLength( minV );
		}
	else
		{
		if (u->getLength() > maxV)
			u->setLength( maxV );
		if (u->getLength() < minV)
			u->setLength( minV );
		}

	/* get downpass sequence if tree topology has changed */
	if (topologyHasChanged == true)
		{
		getDownPassSequence();
		}
	
	//setAllUpdateCls( true );
	//setAllUpdateTis( true );
	//flipAllActiveCls();
	//flipAllActiveTis();
		
	return lnProposalProb;
#	endif

}

#undef DEBUG_TBR
double Tree::updateTbr(void) {
			
#	if defined (DEBUG_TBR)
	cout << "Beginning tree" << endl;
	showNodes(root, 3);
#	endif

	/* parameters that should eventually be set in settings.cpp */
	int reconnectionLimit = 1000;
	double heat = 0.6;
	
	/* set flag to -1 for all nodes */
	for (int i=0; i<numNodes; i++)
		getDownPassNode(i)->setFlag( -1 );

	/* pick a branch at random */
	Node *p;
	do
		{
		p = downPassSequence[(int)(ranPtr->uniformRv()*numNodes)];
		} while ( p->getAnc() == NULL );
	Node *q = p->getAnc();
	if (q->getAnc() == NULL)
		{
		p = q;
		q = p->getLft();
		}
	p->setFlag( 0 );
	q->setFlag( 0 );
	
	/* divide up the tree */
	int numAvailableNodes;
	Node *a, *b, *c, *d, *root1, *root2, *availableNodes[2];
	double oldBrlenConnectionLength1 = 1.0;
	if ( p->getLft() == NULL || p->getRht() == NULL || p->getAnc() == NULL )
		{
		/* we cut at a tip branch */
		if ( p->getLft() == NULL && p->getRht() == NULL )
			{
			/* tip branch pointing up */
			if (q->getLft() == p)
				a = q->getRht();
			else
				a = q->getLft();
			b = q->getAnc();
			if (b->getLft() == q)
				b->setLft(a);
			else
				b->setRht(a);
			a->setAnc(b);
			a->setP( a->getP() + q->getP() );
			a->setFlag( 0 );
			b->setFlag( 1 );
			q->setLft(p);
			q->setRht(NULL);
			q->setAnc(NULL);
			root2 = a;
			while (root2->getAnc() != NULL)
				root2 = root2->getAnc();
			p->setAnc(NULL);
			availableNodes[0] = q;
			numAvailableNodes = 1;
			q->setLft(NULL);
			q->setRht(NULL);
			q->setAnc(NULL);
			q->setP( 0.0 );
			q->setFlag( 0 );
			root1 = p;
			}
		else
			{
			/* tip branch pointing down (at root) */
			root1 = p;
			p->setLft(NULL);
			p->setRht(NULL);
			p->setAnc(NULL);
			p->setFlag( 0 );
			p->setP( q->getP() );
			q->setAnc(NULL);
			a = q->getLft();
			b = q->getRht();
			a->setFlag( 0 );
			b->setFlag( 1 );
			while (b->getRht() != NULL)
				{
				c = b->getRht();
				d = b->getLft();
				
				a->setP( a->getP()+b->getP() );
				b->setP( c->getP() / 2.0 );
				c->setP( c->getP() / 2.0 );

				b->setLft(a);
				b->setRht(d);
				b->setAnc(q);
				a->setAnc(b);
				d->setAnc(b);
				q->setLft(b);
				q->setRht(c);
				c->setAnc(q);
				
				a = q->getLft();
				b = q->getRht();
				}
			b->setLft(a);
			b->setRht(NULL);
			b->setAnc(NULL);
			a->setAnc(b);
			a->setP( a->getP()+b->getP() );
			b->setP( 0.0 );
			root2 = b;
			availableNodes[0] = q;
			numAvailableNodes = 1;
			q->setLft(NULL);
			q->setRht(NULL);
			q->setAnc(NULL);
			q->setP( 0.0 );
			}
		}
	else
		{
		/* cut at an interior branch */
		if (q->getLft() == p)
			a = q->getRht();
		else
			a = q->getLft();
		a->setFlag( 0 );
		b = q->getAnc();
		b->setFlag( 1 );
		if (b->getLft() == q)
			b->setLft(a);
		else
			b->setRht(a);
		a->setAnc(b);
		a->setP( a->getP()+q->getP() );
		root2 = a;
		while (root2->getAnc() != NULL)
			root2 = root2->getAnc();
		availableNodes[0] = q;
		q->setLft(NULL);
		q->setRht(NULL);
		q->setAnc(NULL);
		q->setP( 0.0 );
		p->setAnc(NULL);

		a = p->getLft();
		b = p->getRht();
		a->setFlag( 0 );
		b->setFlag( 1 );
		oldBrlenConnectionLength1 = a->getP() + b->getP();
		while (b->getRht() != NULL)
			{

			c = b->getRht();
			d = b->getLft();
			a->setP( a->getP()+b->getP() );
			b->setP( c->getP() / 2.0 );
			c->setP( c->getP() / 2.0 );

			b->setLft(a);
			b->setRht(d);
			b->setAnc(p);
			a->setAnc(b);
			d->setAnc(b);
			p->setLft(b);
			p->setRht(c);
			c->setAnc(p);
			
			a = p->getLft();
			b = p->getRht();
			}
		b->setLft(a);
		b->setRht(NULL);
		b->setAnc(NULL);
		a->setAnc(b);
		a->setP( a->getP()+b->getP() );
		b->setP( 0.0 );
		root1 = b;
		availableNodes[1] = p;
		numAvailableNodes = 2;
		p->setLft(NULL);
		p->setRht(NULL);
		p->setAnc(NULL);
		p->setP( 0.0 );
		}
	
	/* get down pass sequence for each subtree */
	Node **dp1 = new Node*[numNodes];
	if ( !dp1 )
		{
		std::cerr << "ERROR: Could not allocate dp1" << std::endl;
		exit (0);
		}
	Node **dp2 = new Node*[numNodes];
	if ( !dp2 )
		{
		std::cerr << "ERROR: Could not allocate dp2" << std::endl;
		exit (0);
		}
	int nNodes1 = getDownPassSequence(dp1, root1);
	int nNodes2 = getDownPassSequence(dp2, root2);
	
	/* get the lengths (sum of branch-length proportions) for both subtrees */
	double treeLengthSubTree1 = 0.0, treeLengthSubTree2 = 0.0;
	for (int i=0; i<nNodes1; i++)
		{
		p = dp1[i];
		if (p->getAnc() != NULL)
			treeLengthSubTree1 += p->getP();
		}
	for (int i=0; i<nNodes2; i++)
		{
		p = dp2[i];
		if (p->getAnc() != NULL)
			treeLengthSubTree2 += p->getP();
		}
	
	/* calculate distance from cut portion of the tree (for reconnection limit business) */
	for (int i=0; i<nNodes1; i++)
		{
		p = dp1[i];
		if ( p->getFlag() >= 0 )
			continue;
		if (p->getLft() != NULL)
			{
			if ( p->getLft()->getFlag() >= 0 )
				p->setFlag( p->getLft()->getFlag() + 1 );
			}
		if (p->getRht() != NULL)
			{
			if ( p->getRht()->getFlag() >= 0 )
				p->setFlag( p->getRht()->getFlag() + 1 );
			}
		}
	for (int i=nNodes1-1; i>=0; i--)
		{
		p = dp1[i];
		if (p->getFlag() >= 0)
			continue;
		if (p->getAnc() != NULL)
			p->setFlag( p->getAnc()->getFlag() + 1 );
		if (p->getLft() != NULL)
			{
			if (p->getLft()->getFlag() >= 0)
				p->setFlag( p->getLft()->getFlag() + 1 );
			}
		if (p->getRht() != NULL)
			{
			if (p->getRht()->getFlag() >= 0)
				p->setFlag( p->getRht()->getFlag() + 1 );
			}
		}
	for (int i=0; i<nNodes2; i++)
		{
		p = dp2[i];
		if (p->getFlag() >= 0)
			continue;
		if (p->getLft() != NULL)
			{
			if (p->getLft()->getFlag() >= 0)
				p->setFlag( p->getLft()->getFlag() + 1 );
			}
		if (p->getRht() != NULL)
			{
			if (p->getRht()->getFlag() >= 0)
				p->setFlag( p->getRht()->getFlag() + 1 );
			}
		}
	for (int i=nNodes2-1; i>=0; i--)
		{
		p = dp2[i];
		if (p->getFlag() >= 0)
			continue;
		if (p->getAnc() != NULL)
			p->setFlag( p->getAnc()->getFlag() + 1 );
		if (p->getLft() != NULL)
			{
			if (p->getLft()->getFlag() >= 0)
				p->setFlag( p->getLft()->getFlag() + 1 );
			}
		if (p->getRht() != NULL)
			{
			if (p->getRht()->getFlag() >= 0)
				p->setFlag( p->getRht()->getFlag() + 1 );
			}
		}

#	if defined (DEBUG_TBR)
	cout << "Subtree 1:" << endl;
	showNodes(root1, 3); 
	cout << "Subtree 2:" << endl;
	showNodes(root2, 3); 
#	endif
	
	/* Get the lengths of the branch for the old reconnection point. Note that we cut
	   the tree into two, producing two subtrees rooted at root1 and root2. These
	   trees are largely equivalent. However, the subtree rooted at root2 is the one
	   which holds the branch that was originally the connection point for the
	   subtree rooted at root1; it has a branch that was "healed"---lengthened to take
	   into account the removal of an incident branch. */
	double oldBrlenConnectionLength2 = 0.0;
	for (int i=0; i<nNodes2; i++)
		{
		p = dp2[i];
		if (p->getFlag() == 0)
			{
			oldBrlenConnectionLength2 = p->getP();
			break;
			}
		}
		
	/* initialize state sets on the two subtrees */
	int length1 = modelPtr->getStateSetPtr()->initializeStateSets(&dp1[0], nNodes1);
	int length2 = modelPtr->getStateSetPtr()->initializeStateSets(&dp2[0], nNodes2);
	
	/* allocate matrix to hold parsimony scores for different attachment points */
	int dim1, dim2;
	if (nNodes1 > 1)
		dim1 = nNodes1 - 1;
	else
		dim1 = 1;
	if (nNodes2 > 1)
		dim2 = nNodes2 - 1;
	else
		dim2 = 1;
	int **scoreMatrix = new int*[dim1];
	if (!scoreMatrix)
		{
		std::cerr << "ERROR: Could not allcoate matrix" << std::endl;
		exit (0);
		}
	scoreMatrix[0] = new int[dim1 * dim2];
	if (!scoreMatrix[0])
		{
		std::cerr << "ERROR: Could not allcoate matrix[0]" << std::endl;
		exit (0);
		}
	for (int i=1; i<dim1; i++)
		scoreMatrix[i] = scoreMatrix[i-1] + dim2;		
	for (int i=0; i<dim1; i++)
		for (int j=0; j<dim2; j++)
			scoreMatrix[i][j] = 0;
	double **probs = new double*[dim1];
	if (!probs)
		{
		printf ("   ERROR: Could not allocate probs\n");
		exit (0);
		}
	probs[0] = new double[dim1 * dim2];
	if (!probs[0])
		{
		printf ("   ERROR: Could not allocate probs[0]\n");
		exit (0);
		}
	for (int i=1; i<dim1; i++)
		probs[i] = probs[i-1] + dim2;		
	for (int i=0; i<dim1; i++)
		for (int j=0; j<dim2; j++)
			probs[i][j] = 0.0;
	
	/* calculate parsimony scores for all possible reattachment points */
	double z = (lambda / (4.0/3.0 + lambda));
	double a0 = log(0.25 + 0.75 * z);
	double a1 = log(0.25 - 0.25 * z);
	double marginalLikelihoodTerm = alignmentPtr->getNumChar() * ((2 * numTaxa - 3) * a0 - log(4.0) );
	double tempFactor = a1 - a0;
	double oldLnLike = 0.0;
	for (int i=0; i<nNodes1; i++)
		{
		p = dp1[i];
		if (p->getAnc() != NULL || (p->getAnc() == NULL && nNodes1 == 1))
			{
			for (int j=0; j<nNodes2; j++)
				{
				q = dp2[j];
				if (q->getAnc() != NULL || (q->getAnc() == NULL && nNodes2 == 1))
					{
					if ( p->getFlag() + q->getFlag() < reconnectionLimit )
						{
						unsigned *ss1 = modelPtr->getStateSetPtr()->getStsPtr( 1, p->getIndex() );
						unsigned *ss2 = modelPtr->getStateSetPtr()->getStsPtr( 1, q->getIndex() );
						int length = length1 + length2;
						for (int c=0; c<alignmentPtr->getNumChar(); c++)
							{
							unsigned x = (*ss1);
							unsigned y = (*ss2);
							unsigned zA = (x & y);
							if ( zA == 0 )
								length += alignmentPtr->getNumOfPattern( c );
							ss1++;
							ss2++;
							}
						scoreMatrix[i][j] = length;
						}
					else
						{
						scoreMatrix[i][j] = 1000000000;
						}
					probs[i][j] = tempFactor * scoreMatrix[i][j] + marginalLikelihoodTerm;
					if ( p->getFlag() + q->getFlag() == 0 )
						oldLnLike = probs[i][j];
					}
				}
			}
		}
#	if defined (DEBUG_TBR)
	for (int i=0; i<dim1; i++)
		{
		for (int j=0; j<dim2; j++)
			cout << fixed << setprecision(4) << probs[i][j] << " ";
		cout << endl;
		}
#	endif

	/* rescale log likelihoods, and calculate probabilities for reattachment points */
	double maxLnProb = -100000000000.0;
	for (int i=0; i<dim1; i++)
		{
		for (int j=0; j<dim2; j++)
			{
			probs[i][j] *= heat;
			if (probs[i][j] > maxLnProb)
				maxLnProb = probs[i][j];
			}
		}
	double sum = 0.0;
	for (int i=0; i<dim1; i++)
		{
		for (int j=0; j<dim2; j++)
			{
			probs[i][j] = probs[i][j] - maxLnProb;
			if (probs[i][j] < -150.0)
				probs[i][j] = 0.0;
			else
				probs[i][j] = exp(probs[i][j]);
			sum += probs[i][j];
			}
		}
	for (int i=0; i<dim1; i++)
		{
		for (int j=0; j<dim2; j++)
			{
			probs[i][j] /= sum;
			}
		}
		
	/* pick a reattachment point */
	double ran = ranPtr->uniformRv();
	int i, j;
	for (i=0, sum=0.0; i<dim1; i++)
		{
		for (j=0; j<dim2; j++)
			{
			sum += probs[i][j];
			if (ran < sum)
				break;
			}
		if (ran < sum)
			break;
		}
	p = dp1[i];
	q = dp2[j];
	
	/* get the marginal likelihood of the new reattachment point */
	double newLnLike = tempFactor * scoreMatrix[i][j] + marginalLikelihoodTerm;
	
	/* get the length of the branch for the new reattachment point */
	double newBrlenConnectionLength2 = q->getP();

#	if defined (DEBUG_TBR)
	cout << "i=" << i << " j=" << j << endl;
	cout << "p=" << p->getIndex() << " q=" << q->getIndex() << endl;
	cout << "length1 = " << length1 << endl;
	cout << "length2 = " << length2 << endl;
	cout << "dp1 = " << dp1 << endl;
	cout << "dp2 = " << dp2 << endl;
	cout << "scoreMatrix = " << scoreMatrix << endl;
	cout << "scoreMatrix[0] = " << scoreMatrix[0] << endl;
	cout << "probs = " << probs << endl;
	cout << "probs[0] = " << probs[0] << endl;
	cout << "this = " << this << endl;
	
	for (int i=0; i<dim1; i++)
		{
		for (int j=0; j<dim2; j++)
			cout << setw(4) << scoreMatrix[i][j] << " ";
		cout << endl;
		}
	
	for (int i=0; i<dim1; i++)
		{
		for (int j=0; j<dim2; j++)
			cout << fixed << setprecision(4) << probs[i][j] << " ";
		cout << endl;
		}
#	endif
		
	/* link up the two subtrees at branches p & q */
	double newBrlenConnectionLength1 = 1.0;
	if (numAvailableNodes == 2)
		{
		/* reroot subtree 1 such that p is to the left or right of the root (node a is the root; tree is rooted) */
		newBrlenConnectionLength1 = p->getP();
		a = availableNodes[0];
		b = root1;
		c = b->getLft();
		for (int i=0; i<nNodes1; i++)
			dp1[i]->setMarked(false);
		markBranchesDown(p);
		a->setLft(c);
		a->setRht(b);
		b->setAnc(a);
		c->setAnc(a);
		b->setMarked(false);
		b->setLft(NULL);
		b->setRht(NULL);
		a->setAnc(NULL);
		double u = ranPtr->uniformRv() * c->getP();
		b->setP( c->getP()-u );
		c->setP( u );
		a->setP( 0.0 );
		
		while ( (a->getLft()->getMarked() == true && a->getLft() != p) || (a->getRht()->getMarked() == true && a->getRht() != p) )
			{
			Node *e;
			if ( a->getLft()->getMarked() == true && a->getRht()->getMarked() == true )
				{
				std::cerr << "ERROR: Marked to the left and the right" << std::endl;
				exit(1);
				}
			if (a->getLft()->getMarked() == true)
				{
				b = a->getLft();
				c = a->getRht();
				}
			else
				{
				b = a->getRht();
				c = a->getLft();
				}
			if (b->getLft() == NULL && b->getRht() == NULL)
				{
				if (b != p)
					{
					std::cerr << "ERROR: p is not the last marked node" << std::endl;
					exit(1);
					}
				break;
				}
			else
				{
				if (b->getLft()->getMarked() == true)
					{
					d = b->getLft();
					e = b->getRht();
					}
				else
					{
					d = b->getRht();
					e = b->getLft();
					}
				}
			b->setLft(e);
			b->setRht(c);
			e->setAnc(b);
			c->setAnc(b);
			a->setLft(b);
			a->setRht(d);
			b->setAnc(a);
			d->setAnc(a);
			b->setMarked(false);
			c->setP( c->getP()+b->getP() );
			u = ranPtr->uniformRv() * d->getP();
			b->setP( d->getP()-u );
			d->setP( u );
			}
		b = availableNodes[1];
		c = q->getAnc();
		u = ranPtr->uniformRv() * q->getP();
		b->setP( u );
		q->setP( q->getP()-u );
		
		if (c->getLft() == q)
			{
			c->setLft(b);
			b->setAnc(c);
			b->setLft(q);
			q->setAnc(b);
			b->setRht(a);
			a->setAnc(b);
			}
		else
			{
			c->setRht(b);
			b->setAnc(c);
			b->setRht(q);
			q->setAnc(b);
			b->setLft(a);
			a->setAnc(b);
			}
		a->setP( 1.0 - treeLengthSubTree1 - treeLengthSubTree2 );
		root = root2;
		}
	else
		{
		a = availableNodes[0];
		if (nNodes1 == 1)
			{
			c = q->getAnc();
			if (c->getLft() == q)
				{
				c->setLft(a);
				a->setAnc(c);
				a->setLft(q);
				q->setAnc(a);
				a->setRht(p);
				p->setAnc(a);
				}
			else
				{
				c->setRht(a);
				a->setAnc(c);
				a->setRht(q);
				q->setAnc(a);
				a->setLft(p);
				p->setAnc(a);
				}
			double u = ranPtr->uniformRv() * q->getP();
			a->setP( u );
			q->setP( q->getP()-u );
			root = root2;
			}
		else if (nNodes2 == 1)
			{
			c = p->getAnc();
			if (c->getLft() == p)
				{
				c->setLft(a);
				a->setAnc(c);
				a->setLft(p);
				p->setAnc(a);
				a->setRht(q);
				q->setAnc(a);
				}
			else
				{
				c->setRht(a);
				a->setAnc(c);
				a->setRht(p);
				p->setAnc(a);
				a->setLft(q);
				q->setAnc(a);
				}
			double u = ranPtr->uniformRv() * p->getP();
			a->setP( u );
			p->setP( p->getP()-u );
			root = root1;
			}
		else
			{
			std::cerr << "ERROR: Why am I here?" << std::endl;
			exit(1);
			}
		}

	/* update downpass sequence */
	getDownPassSequence();
	
	/* calculate proposal ratio */
	double lnProposalProb = heat * (oldLnLike - newLnLike) + 
	                        (log(newBrlenConnectionLength2) - log(oldBrlenConnectionLength2)) + 
							(log(newBrlenConnectionLength1) - log(oldBrlenConnectionLength1));
	
	/* check all of the branch lengths, to make certain that they are a minimum length */
	bool resetBrProp = false;
	sum = 0.0;
	for (int n=0; n<numNodes; n++)
		{
		Node *p = getDownPassNode(n);
		if (p->getAnc() != NULL)
			{
			if ( p->getP() < BRPROP_MIN )
				{
				p->setP( BRPROP_MIN );
				resetBrProp = true;
				}
			sum += p->getP();
			}
		}
	if ( resetBrProp == true )
		{
		for (int n=0; n<numNodes; n++)
			{
			Node *p = getDownPassNode(n);
			if (p->getAnc() != NULL)
				p->setP( p->getP() / sum );
			}
		}
	
#	if defined (DEBUG_TBR)
	for (int n=0; n<numNodes; n++)
		cout << setw(3) << downPassSequence[n]->getIndex();
	cout << endl;
	cout << "oldLnLike = " << oldLnLike << endl;
	cout << "newLnLike = " << newLnLike << endl;
	cout << "oldBrlenConnectionLength1 = " << oldBrlenConnectionLength1 << endl;
	cout << "newBrlenConnectionLength1 = " << newBrlenConnectionLength1 << endl;
	cout << "oldBrlenConnectionLength2 = " << oldBrlenConnectionLength2 << endl;
	cout << "newBrlenConnectionLength2 = " << newBrlenConnectionLength2 << endl;
	cout << "Tree length = " << treeLength() << endl;
	cout << "Ending tree" << endl;
	showNodes(root, 3);
#	endif

	/* free memory */
	delete [] dp1;
	delete [] dp2;
	delete [] scoreMatrix[0];
	delete [] scoreMatrix;
	delete [] probs[0];
	delete [] probs;

	return lnProposalProb;
}

void Tree::clone(Tree &t) {

	if ( numNodes != t.numNodes )
		{
		if (nodes != NULL)
			{
			delete [] nodes;
			}
		if (downPassSequence != NULL)
			{
			delete [] downPassSequence;
			}
		nodes = new Node[t.numNodes];
		for (int i=0; i<t.numNodes; i++)
			nodes[i].setIndex(i);
		downPassSequence = new Node*[t.numNodes];
		for (int i=0; i<t.numNodes; i++)
			downPassSequence[i] = NULL;
		}
		
	lambda         = t.lambda;
	alpha0         = t.alpha0;
	numTaxa        = t.numTaxa;
	numNodes       = t.numNodes;
	alignmentPtr   = t.alignmentPtr;
	isTreeFixed    = t.isTreeFixed;

	for (int i=0; i<numNodes; i++)
		{
		Node *p = &t.nodes[i];
		int idx = p->getIndex();
		Node *q = &nodes[idx];
		
		q->setP( p->getP() );
		q->setIndex( p->getIndex() );
		q->setName( p->getName() );
		q->setIsLeaf( p->getIsLeaf() );
		
		Node *a = p->getLft();
		if (a != NULL)
			{
			idx = a->getIndex();
			Node *b = &nodes[idx];
			q->setLft(b);
			}
		else
			q->setLft(NULL);

		a = p->getRht();
		if (a != NULL)
			{
			idx = a->getIndex();
			Node *b = &nodes[idx];
			q->setRht(b);
			}
		else
			q->setRht(NULL);
			
		a = p->getAnc();
		if (a != NULL)
			{
			idx = a->getIndex();
			Node *b = &nodes[idx];
			q->setAnc(b);
			}
		else
			q->setAnc(NULL);
		}
	int idx = t.root->getIndex();
	root = &nodes[idx];
	for (int i=0; i<numNodes; i++)
		{
		Node *p = t.downPassSequence[i];
		idx = p->getIndex();
		Node *q = &nodes[idx];
		downPassSequence[i] = q;
		}

}

int Tree::dex(Node *p) {

	if ( p == NULL )
		return -1;
	else
		return p->getIndex();

}

void Tree::getDownPassSequence(void) {

	int i = 0;
	passDown(root, &i);
	
}

int Tree::getDownPassSequence(Node **dp, Node *r) {

	int i = 0;
	passDown(r, dp, &i);
	return i;
}

void Tree::passDown(Node *p, Node **dp, int *x) {

	if (p != NULL)
		{
		passDown(p->getLft(), dp, x);
		passDown(p->getRht(), dp, x);
		dp[(*x)++] = p;
		}
}

void Tree::markBranchesDown(Node *p) {

	Node *r = p;
	while (r->getAnc() != NULL)
		{
		r->setMarked(true);
		r = r->getAnc();
		}
}

void Tree::passDown(Node *p, int *x) {

	if (p != NULL)
		{
		passDown(p->getLft(), x);
		passDown(p->getRht(), x);
		downPassSequence[(*x)++] = p;
		}
		
}

void Tree::print(void) {

	showNodes(root, 3);
	
}

void Tree::showNodes(Node *p, int indent) {

	if (p != NULL)
		{
		for (int i=0; i<indent; i++)
			std::cout << " ";
		std::cout << dex(p) << " (" << dex(p->getLft()) << ", " << dex(p->getRht()) << ", " << dex(p->getAnc()) << ") " << std::fixed << std::setprecision(5) << p->getP();
		if (p->getIsLeaf() == true )
			std::cout << " (" << p->getName() << ") ";
		if (p == root)
			std::cout << " <- Root" << std::endl;
		else
			std::cout << std::endl;
		showNodes (p->getLft(),  indent + 2);
		showNodes (p->getRht(), indent + 2);
		}
   
}
