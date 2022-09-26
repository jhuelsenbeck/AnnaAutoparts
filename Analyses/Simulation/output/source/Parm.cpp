#include "Model.h"
#include "Parm.h"
#include "ParmAsrv.h"
#include "ParmFreqs.h"
#include "ParmLength.h"
#include "ParmSubrates.h"
#include "ParmTree.h"
#include <iostream>



Parm::Parm(MbRandom *rp, Model *mp, std::string nm)  {

	ranPtr = rp;
	parmName = nm;
	modelPtr = mp;
	
}

Parm::~Parm(void) {

}

Parm &Parm::operator=(Parm &b) {

	if (this != &b)
		{
		/* copy base class data */
		ranPtr   = b.ranPtr;
		parmName = b.parmName;
		modelPtr = b.modelPtr;
		
		/* We need to downcast the object to the derived class pointers.
		   This allows us to call the correct assignment functions in the
		   derived class. */
		
		/* check to see if the object is a substitution rate parameter */
		{
		SubRates *thatDerivedPtr = dynamic_cast<SubRates *>(&b);
		SubRates *thisDerivedPtr = dynamic_cast<SubRates *>(this);
		if ( thatDerivedPtr != 0 && thisDerivedPtr != 0 )
			{
			thisDerivedPtr->clone( *thatDerivedPtr );
			goto exitOperator;
			}
		}

		/* check to see if the object is a gamma shape parameter */
		{
		Asrv *thatDerivedPtr = dynamic_cast<Asrv *>(&b);
		Asrv *thisDerivedPtr = dynamic_cast<Asrv *>(this);
		if ( thatDerivedPtr != 0 && thisDerivedPtr != 0 )
			{
			thisDerivedPtr->clone( *thatDerivedPtr );
			goto exitOperator;
			}
		}

		/* check to see if the object is a base frequency parameter */
		{
		BaseFreqs *thatDerivedPtr = dynamic_cast<BaseFreqs *>(&b);
		BaseFreqs *thisDerivedPtr = dynamic_cast<BaseFreqs *>(this);
		if ( thatDerivedPtr != 0 && thisDerivedPtr != 0 )
			{
			thisDerivedPtr->clone( *thatDerivedPtr );
			goto exitOperator;
			}
		}

		/* check to see if the object is a tree parameter */
		{
		Tree *thatDerivedPtr = dynamic_cast<Tree *>(&b);
		Tree *thisDerivedPtr = dynamic_cast<Tree *>(this);
		if ( thatDerivedPtr != 0 && thisDerivedPtr != 0 )
			{
			thisDerivedPtr->clone( *thatDerivedPtr );
			goto exitOperator;
			}
		}

		/* check to see if the object is a tree parameter */
		{
		TreeLength *thatDerivedPtr = dynamic_cast<TreeLength *>(&b);
		TreeLength *thisDerivedPtr = dynamic_cast<TreeLength *>(this);
		if ( thatDerivedPtr != 0 && thisDerivedPtr != 0 )
			{
			thisDerivedPtr->clone( *thatDerivedPtr );
			goto exitOperator;
			}
		}
		
		std::cout << "Problem in Parameter assignment operator" << std::endl;
		exit(1);
			
		exitOperator:
			;
		}
	return *this;

}
