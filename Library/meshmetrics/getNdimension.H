#include "mesh.h"

dint getNdimension(dint fieldType);

dint getNdimension(dint fieldType){

	dint n = 0;

	switch (fieldType)
	{
	case volumeScalar: 		n = 1; break;
	case volumeScalar2: 	n = 1; break;
	case volumeScalar3: 	n = 1; break;
	case surfaceScalar: 	n = 1; break;
	case vertexScalar: 		n = 1; break;
	case volumeVector:  	n = 3; break;
	case volumeVector2: 	n = 3; break;
	case volumeVector3: 	n = 3; break;
	case surfaceVector: 	n = 3; break;
	case vertexVector:  	n = 3; break;
	case scalar:    		n = 1; break;
	case vector:    		n = 3; break;
	}

	return n;
}

