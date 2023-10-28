#include "mesh.h"

dfloat dotProduct(dfloat *A,dfloat *B);

dfloat dotProduct(dfloat *A,dfloat *B){

	// returns square of a number

	dfloat AdotB = A[0]*B[0]+A[1]*B[1]+A[2]*B[2];

	return AdotB;

}
