#include "mesh.h"

dfloat magVector(dfloat *X);

dfloat magVector(dfloat *X){

	dfloat magnitudeVector = 0.0;

	dfloat a = X[0];
	dfloat b = X[1];
	dfloat c = X[2];

	magnitudeVector = sqrt(a*a +b*b +c*c);

	return magnitudeVector;
}
