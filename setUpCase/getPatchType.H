#include "mesh.h"

const char* getPatchType(dint variableName);

const char* getPatchType(dint variableName){

	switch (variableName)
	{
	case periodicPatch:   return "periodicPatch";
	case boundaryPatch:   return "boundaryPatch";
	case emptyPatch:      return "emptyPatch";
	}

	return "none";
}
