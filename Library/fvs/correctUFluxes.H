#include "mesh.h"

//#include "../Library/fvs/correctUFluxes.H"

void correctUFluxes(mesh_t *mesh);

void correctUVolumetricFluxes(mesh_t *mesh);
void correctUMassFluxes(mesh_t *mesh);


void correctUFluxes(mesh_t *mesh){

	correctUVolumetricFluxes(mesh);

	correctUMassFluxes(mesh);

}


void correctUVolumetricFluxes(mesh_t *mesh){

	for (int cell = 0; cell < mesh->nicells; cell++) {

		for (int face = 0; face < mesh->ncellFaces; face++) {

			dint fid 	= face + cell*mesh->ncellFaces;
			int cell_n	=  mesh->faceToNeighborCell[fid];

			mesh->phiU[fid] = mesh->phiHbyA[fid] - mesh->phigradPbyA[fid];

		}

	}



}



void correctUMassFluxes(mesh_t *mesh){


	for( dint cell = 0; cell < mesh->ncells;cell++){
		for (dint face = 0; face < mesh->ncellFaces; ++face) {
			dint fid = face + cell*mesh->ncellFaces;
			mesh->rhophiU[fid] = mesh->phiU[fid]*mesh->rhof[fid];
		}
	}


}
