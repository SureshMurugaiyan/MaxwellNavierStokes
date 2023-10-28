#include "mesh.h"

void generateComputationalDomain(mesh_t *mesh);

void generateComputationalDomain(mesh_t *mesh){


	for (dint i = 0; i < mesh->nverticesx; i++){

		dfloat minVal = 0.0;
		dfloat maxVal = 1.0;

		dfloat delta=((( maxVal - minVal)/(mesh->nverticesx-1)));

		mesh->xi[i] =  minVal + delta*i;

	}


	for (dint j = 0; j < mesh->nverticesy; j++){

		dfloat minVal = 0.0;
		dfloat maxVal = 1.0;

		dfloat delta=((( maxVal - minVal)/(mesh->nverticesy-1)));

		mesh->eta[j] =  minVal + delta*j;

	}


	for (dint k = 0; k < mesh->nverticesz; k++){

		dfloat minVal = 0.0;
		dfloat maxVal = 1.0;

		dfloat delta=((( maxVal - minVal)/(mesh->nverticesz-1)));

		mesh->zeta[k] =  minVal + delta*k;

	}


}
