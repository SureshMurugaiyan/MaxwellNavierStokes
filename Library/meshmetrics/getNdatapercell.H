#include "mesh.h"

dint getNdatapercell(dint fieldType);

dint getNdatapercell(dint fieldType){

	dint n = 0;

	switch (fieldType)
	{
	case volumeScalar: 	n = 1; break;
	case volumeScalar2: n = 1; break;
	case volumeScalar3: n = 1; break;
	case surfaceScalar: n = 4; break;
	case vertexScalar: 	n = 4; break;
	case volumeVector:  n = 1; break;
	case volumeVector2: n = 1; break;
	case volumeVector3: n = 1; break;
	case surfaceVector: n = 4; break;
	case vertexVector:  n = 4; break;
	case scalar:    n = 1; break;
	case vector:    n = 1; break;
	}

	return n;
}
