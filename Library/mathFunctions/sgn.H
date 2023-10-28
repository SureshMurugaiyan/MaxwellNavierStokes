#include "mesh.h"

dfloat sgn(dfloat F);

dfloat sgn(dfloat F){

	// returns 1 if Flux is zero or greater than zero

	int sgn = 1;

	if (F < (-1.0e-15)){
		sgn =  -1;
	}

	return sgn;

}

