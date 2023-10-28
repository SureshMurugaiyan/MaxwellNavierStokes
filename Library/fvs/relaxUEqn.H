#include "mesh.h"

//#include "../Library/fvs/relaxUEqn.H"

void relaxUEqn(mesh_t *mesh);

void relaxUEqn(mesh_t *mesh){



	if (mesh->currentNCorrector == 1) {

		return;
	}



	for (int cell = 0; cell < mesh->nicells; cell++) {
		for (dint dimension = 0; dimension < mesh->ndimension; dimension++) {

			dint cell_dim = dimension + cell*mesh->ndimension;

			mesh->BP[cell_dim]  =  mesh->BP[cell_dim] + (1-mesh->UEqnlambda)*(mesh->AP[cell_dim]/mesh->UEqnlambda)*mesh->Ustar[cell_dim];

		}
	}


	for (int cell = 0; cell < mesh->nicells; cell++) {
		for (dint dimension = 0; dimension < mesh->ndimension; dimension++) {

			dint cell_dim = dimension + cell*mesh->ndimension;

			mesh->AP[cell_dim]  =  mesh->AP[cell_dim] /mesh->UEqnlambda ;

		}
	}

}
