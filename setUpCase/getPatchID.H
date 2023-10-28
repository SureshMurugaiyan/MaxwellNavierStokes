#include "mesh.h"
#include "BoundaryNames.H"

const char* getPatchID(dint variableName);

const char* getPatchID(dint variableName){

	switch (variableName)
	{
	case topBoundary: 	 return "topBoundary";
	case bottomBoundary: return "bottomBoundary";
	case leftBoundary: 	 return "leftBoundary";
	case rightBoundary:  return "rightBoundary";
	case frontBoundary:  return "frontBoundary";
	case backBoundary:   return "backBoundary";
	}

	return "none";
}
