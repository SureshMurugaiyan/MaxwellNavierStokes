#include "mesh.h"

dint computeLinearIndex(dint i, dint j, dint k, dint nx, dint ny, dint nz);

dint computeLinearIndex(dint i, dint j, dint k, dint nx, dint ny, dint nz){

	dint linearIndex =  k + j*nz + i*ny*nz;

	return  linearIndex;


}


