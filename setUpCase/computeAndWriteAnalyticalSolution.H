#include "mesh.h"

//#include "../setUpCase/computeAndWriteAnalyticalSolution.H"

void computeAndWriteAnalyticalSolution(mesh_t *mesh);

void computeAndWriteAnalyticalSolution(mesh_t *mesh){


	/*******************************************************************************************
	 * Template to Code the exact solution and print at each time level
	 * ****************************************************************************************/

	dfloat *phiScalarExact;
	dfloat *phiVectorExact;

	createField(&phiScalarExact,volumeScalar,mesh);
	createField(&phiVectorExact,volumeVector,mesh);

	dfloat t = mesh->time;

	for(int cell = 0;cell < mesh->ncells; ++cell) {

		dint idx	= mesh->dim_x + cell*mesh->ndimension;
		dint idy	= mesh->dim_y + cell*mesh->ndimension;
		dint idz	= mesh->dim_z + cell*mesh->ndimension;

		dfloat x = mesh->C[idx];
		dfloat y = mesh->C[idy];
		dfloat z = mesh->C[idz];

		phiScalarExact[cell] = 0.0;

		phiVectorExact[idx]	 = 0.0;
		phiVectorExact[idy]	 = 0.0;
		phiVectorExact[idz]	 = 0.0;

	}


	// write exact solution

	strcpy(mesh->variableName, "phiScalarExact");
	writeResultsToFile(phiScalarExact,volumeScalar, mesh->L_0,1.0,mesh);

	strcpy(mesh->variableName, "phiVectorExact");
	writeResultsToFile(phiVectorExact,volumeVector, mesh->L_0,1.0,mesh);

	free(phiScalarExact);
	free(phiVectorExact);


}

