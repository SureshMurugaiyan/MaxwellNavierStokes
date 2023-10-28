#include "mesh.h"
#include "BoundaryNames.H"
#include "FieldNames.H"

//#include "../setUpCase/getBoundaryConditionType.H"

dint getBoundaryConditionType(dint patch, dint phiName, mesh_t *mesh);

dint getBoundaryConditionType(dint patch, dint phiName, mesh_t *mesh){

	dint patchID = mesh->patchID[patch];
	dint bctype  = none;

	switch (patchID)
	{
	case topBoundary:
		switch (phiName)
		{
		case U:
			bctype			=	fixedValue;
			break;
		case p:
			bctype			=	zeroNormalGradient;
			break;
		case alpha:
			bctype			=	zeroNormalGradient;
			break;
		case phi:
			bctype			=	zeroNormalGradient;
			break;
		case A:
			bctype			=	fixedValue;
			break;
		}
		break;
		/*******************************************************/
	case bottomBoundary:
		switch (phiName)
		{
		case U:
			bctype			=	fixedValue;
			break;
		case p:
			bctype			=	zeroNormalGradient;
			break;
		case alpha:
			bctype			=	zeroNormalGradient;
			break;
		case phi:
			bctype			=	zeroNormalGradient;
			break;
		case A:
			bctype			=	fixedValue;
			break;
		}
		break;
		/*******************************************************/
	case frontBoundary:
		switch (phiName)
		{
		case U:
			bctype			=	fixedValue;
			break;
		case p:
			bctype			=	zeroNormalGradient;
			break;
		case alpha:
			bctype			=	zeroNormalGradient;
			break;
		case phi:
			bctype			=	fixedValue;
			break;
		case A:
			bctype			=	fixedValue;
			break;
		}
		break;
		/*******************************************************/
	case backBoundary:
		switch (phiName)
		{
		case U:
			bctype			=	fixedValue;
			break;
		case p:
			bctype			=	zeroNormalGradient;
			break;
		case alpha:
			bctype			=	zeroNormalGradient;
			break;
		case phi:
			bctype			=	fixedValue;
			break;
		case A:
			bctype			=	fixedValue;
			break;
		}
		break;
		/*******************************************************/
	case leftBoundary:
			switch (phiName)
			{
			case U:
				bctype			=	zeroNormalGradient;
				break;
			case p:
				bctype			=	zeroNormalGradient;
				break;
			case alpha:
				bctype			=	zeroNormalGradient;
				break;
			case phi:
				bctype			=	zeroNormalGradient;
				break;
			case A:
				bctype			=	zeroNormalGradient;
				break;
			}
			break;
			/*******************************************************/
	 case rightBoundary:
			switch (phiName)
			{
			case U:
				bctype			=	zeroNormalGradient;
				break;
			case p:
				bctype			=	zeroNormalGradient;
				break;
			case alpha:
				bctype			=	zeroNormalGradient;
				break;
			case phi:
				bctype			=	zeroNormalGradient;
				break;
			case A:
				bctype			=	zeroNormalGradient;
				break;
			}
			break;
			/*******************************************************/
	}


	if ( (bctype  == periodic) &&
		 (mesh->periodicNeighborPatchID[patch] == none))
	{
		printf ( "\n" );
		printf ( "Fatal error!\n" );
		printf ( "Periodic Neighbor not set for %15s \n",getPatchID(patchID));
		exit ( 1 );
	}


//	printf ( " %5s: %20s %20s \n",getVariableID(phiName),getPatchID(patchID), getBctypeID(bctype) );


	return  bctype;
}
