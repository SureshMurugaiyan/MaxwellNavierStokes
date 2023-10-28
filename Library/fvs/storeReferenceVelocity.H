#include "mesh.h"

//#include "../Library/fvs/storeReferenceVelocity.H"

void storeReferenceVelocity(mesh_t *mesh);

void storeReferenceVelocity(mesh_t *mesh){

	if (mesh->currentNCorrector==1) {

		for (dint cell = 0; cell < mesh->ncells; cell++) {

			dint cell_O  = cell + mesh->ncells;
			dint cell_OO = cell + mesh->ncells*2;

			for (dint dimension = 0; dimension < mesh->ndimension; dimension++) {

				dint cell_dim 		= dimension + cell*mesh->ndimension;
				dint cell_O_dim 	= dimension + cell_O*mesh->ndimension;
				dint cell_OO_dim 	= dimension + cell_OO*mesh->ndimension;

				mesh->Ustar[cell_dim]   = 0.0;

			}
		}
	}

	if (mesh->currentNCorrector > 1) {

		for (dint cell = 0; cell < mesh->ncells; cell++) {

			dint cell_O  = cell + mesh->ncells;
			dint cell_OO = cell + mesh->ncells*2;

			for (dint dimension = 0; dimension < mesh->ndimension; dimension++) {

				dint cell_dim 		= dimension + cell*mesh->ndimension;
				dint cell_O_dim 	= dimension + cell_O*mesh->ndimension;
				dint cell_OO_dim 	= dimension + cell_OO*mesh->ndimension;

				mesh->Ustar[cell_dim]   = mesh->U[cell_dim];

			}
		}

	}


//	strcpy(mesh->variableName, "Ustar");
//	printVolumeVector(mesh->Ustar, mesh);

}
