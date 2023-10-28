#include "mesh.h"

//#include "../setUpCase/setUpExternalMagneticField.H"

void setUpExternalMagneticField(mesh_t *mesh);

void setUpExternalMagneticField(mesh_t *mesh){

	// Fill the cell center values
	for(int cell = 0; cell < mesh->ncells; ++cell){

		dint idx	= mesh->dim_x + cell*mesh->ndimension;
		dint idy	= mesh->dim_y + cell*mesh->ndimension;
		dint idz	= mesh->dim_z + cell*mesh->ndimension;

		dfloat x = mesh->C[idx];
		dfloat y = mesh->C[idy];
		dfloat z = mesh->C[idz];

		dfloat B0     = 1.0;
		dfloat alphaB = 1.0;

		mesh->B[idx]		=  0.0;

		mesh->B[idy]		=  B0*(1+ alphaB *x);

		mesh->B[idz]		=   0.0;

	}


//	strcpy(mesh->variableName, "B");
//	printVolumeVector(mesh->B, mesh);

}

