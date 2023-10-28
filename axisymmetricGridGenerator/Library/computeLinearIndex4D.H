#include "mesh.h"

dint computeLinearIndex4D(linearIndex4D_t *linearIndex4D);

dint computeLinearIndex4D(linearIndex4D_t *linearIndex4D){

	dint i = linearIndex4D->i;
	dint j = linearIndex4D->j;
	dint k = linearIndex4D->k;
	dint l = linearIndex4D->l;

	dint iN = linearIndex4D->iN;
	dint jN = linearIndex4D->jN;
	dint kN = linearIndex4D->kN;
	dint lN = linearIndex4D->lN;

	dint linearIndex = l + k*lN + j*lN*kN + i*lN*kN*jN;

	return  linearIndex;

}


