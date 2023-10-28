#include "mesh.h"
//#include "../setUpCase/getFieldType.H"

#include "FieldNames.H"

dint getFieldType(dint phiName);

dint getFieldType(dint phiName){

	dint n = volumeScalar;

	switch (phiName)
	{
	case U:			n = volumeVector3;  break;
	case p: 		n = volumeScalar3;  break;
	case alpha: 	n = volumeScalar3;  break;
	}

	return n;
}
