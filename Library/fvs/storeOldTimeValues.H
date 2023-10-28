#include "mesh.h"
//#include "../Library/fvs/storeOldTimeValues.H"

void storeOldTimeValues(mesh_t *mesh);

void storeOldTimeValuesFirstTimeStep(mesh_t *mesh);
void storeOldTimeValuesAfterFirstTimeStep(mesh_t *mesh);


void storeOldTimeValues(mesh_t *mesh){

	if (mesh->timeStep==0) {

		storeOldTimeValuesFirstTimeStep(mesh);

	} else{

		storeOldTimeValuesAfterFirstTimeStep(mesh);

	}


}


void storeOldTimeValuesFirstTimeStep(mesh_t *mesh){

	for (dint cell = 0; cell < mesh->ncells; cell++) {

		dint cell_O  = cell + mesh->ncells;
		dint cell_OO = cell + mesh->ncells*2;

		for (dint dimension = 0; dimension < mesh->ndimension; dimension++) {

			dint cell_dim 		= dimension + cell*mesh->ndimension;
			dint cell_O_dim 	= dimension + cell_O*mesh->ndimension;
			dint cell_OO_dim 	= dimension + cell_OO*mesh->ndimension;

			mesh->U[cell_O_dim]  = mesh->U[cell_dim];
			mesh->U[cell_OO_dim] = mesh->U[cell_dim];

		}

		mesh->p[cell_O]   = mesh->p[cell];
		mesh->p[cell_OO]  = mesh->p[cell];

		mesh->alpha[cell_O]  = mesh->alpha[cell];
		mesh->alpha[cell_OO] = mesh->alpha[cell];

	}

}

void storeOldTimeValuesAfterFirstTimeStep(mesh_t *mesh){

	for (dint cell = 0; cell < mesh->ncells; cell++) {

		dint cell_O  = cell + mesh->ncells;
		dint cell_OO = cell + mesh->ncells*2;

		for (dint dimension = 0; dimension < mesh->ndimension; dimension++) {

			dint cell_dim 		= dimension + cell*mesh->ndimension;
			dint cell_O_dim 	= dimension + cell_O*mesh->ndimension;
			dint cell_OO_dim 	= dimension + cell_OO*mesh->ndimension;

			mesh->U[cell_OO_dim] = mesh->U[cell_O_dim];
			mesh->U[cell_O_dim]  = mesh->U[cell_dim];

		}

		mesh->p[cell_OO] = mesh->p[cell_O];
		mesh->p[cell_O]  = mesh->p[cell];

		mesh->alpha[cell_OO] = mesh->alpha[cell_O];
		mesh->alpha[cell_O]  = mesh->alpha[cell];

	}

}
