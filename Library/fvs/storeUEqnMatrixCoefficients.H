#include "mesh.h"

//#include "../Library/fvs/storeUEqnMatrixCoefficients.H"

void storeUEqnMatrixCoefficients(mesh_t *mesh);

void storeUEqnMatrixCoefficients(mesh_t *mesh){

	for (int cell = 0; cell < mesh->nicells; cell++) {
		for (dint dimension = 0; dimension < mesh->ndimension; dimension++) {

			dint cell_dim = dimension + cell*mesh->ndimension;

			mesh->AAP[cell_dim]  =  mesh->AP[cell_dim];
			mesh->BBP[cell_dim]  =  mesh->BP[cell_dim];

		}
	}


	for (int cell = 0; cell < mesh->nicells; cell++) {
		for (int face = 0; face < mesh->ncellFaces; face++) {
			for (int dim = 0; dim < mesh->ndimension; dim++) {

				dint fid 	 = face + cell*mesh->ncellFaces;
				dint fid_dim = dim + fid*mesh->ndimension;

				mesh->AAN[fid_dim]  =  mesh->AN[fid_dim];

			}
		}
	}



}
