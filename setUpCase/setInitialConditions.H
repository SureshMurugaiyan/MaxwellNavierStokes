#include "mesh.h"
#include "FieldNames.H"
#include "DomainNames.H"

//#include "../setUpCase/setInitialConditions.H"

void setInitialConditions(mesh_t *mesh);

void setInitialConditions(mesh_t *mesh){

	dfloat t = mesh->time;

	// Fill the cell center values
	for(int cell = 0; cell < mesh->ncells; ++cell){

		dint idx	= mesh->dim_x + cell*mesh->ndimension;
		dint idy	= mesh->dim_y + cell*mesh->ndimension;
		dint idz	= mesh->dim_z + cell*mesh->ndimension;

		dfloat x = mesh->C[idx];
		dfloat y = mesh->C[idy];
		dfloat z = mesh->C[idz];

		mesh->U[idx]		=  0.0; // x-velocity
		mesh->U[idy]		=  0.0; // y-velocity
		mesh->U[idz]		=  0.0; // z-velocity
		mesh->p[cell]		=  1.0;


		if(mesh->domain[cell] == Gas){
			mesh->alpha[cell]	=  mesh->alphaG;
			mesh->rho[cell]		=  mesh->rhoG;
			mesh->mu[cell]		=  mesh->muG;
			mesh->sigma[cell]	=  mesh->sigmaG;

		}


		if(mesh->domain[cell] == Liquid){
			mesh->alpha[cell]	=  mesh->alphaL;
			mesh->rho[cell]		=  mesh->rhoL;
			mesh->mu[cell]		=  mesh->muL;
			mesh->sigma[cell]	=  mesh->sigmaL;
		}

	}


//	strcpy(mesh->variableName, "U");
//	printVolumeVector(mesh->U, mesh);
//
//	strcpy(mesh->variableName, "p");
//	printVolumeScalar(mesh->p, mesh);

//	strcpy(mesh->variableName, "alpha");
//	printVolumeScalar(mesh->alpha, mesh);

//	strcpy(mesh->variableName, "rho");
//	printVolumeScalar(mesh->rho, mesh);

//	strcpy(mesh->variableName, "mu");
//	printVolumeScalar(mesh->mu, mesh);

}

