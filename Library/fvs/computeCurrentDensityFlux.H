#include "mesh.h"

//#include "../Library/fvs/computeCurrentDensityFlux.H"

void computeCurrentDensityFlux(mesh_t *mesh);

void computeCurrentDensityFlux(mesh_t *mesh){



	for (int cell = 0; cell < mesh->nicells; cell++) {
		for (int face = 0; face < mesh->ncellFaces; face++) {

			dint fid 	= face + cell*mesh->ncellFaces;

			dint OCell 	= mesh->faceToOwnerCell[fid];
			dint NCell 	= mesh->faceToNeighborCell[fid];

			dfloat Flux = - (mesh->phi[NCell]-mesh->phi[OCell])*mesh->rdelta[fid]*mesh->magSf[fid];

			mesh->sigmagradPhiFlux[fid] =	mesh->sigmaf[fid]*Flux;

		}
	}


	for (int cell = 0; cell < mesh->nicells; cell++) {
		for (int face = 0; face < mesh->ncellFaces; face++) {

			dint fid 	= face + cell*mesh->ncellFaces;

			mesh->phiJ[fid] =	mesh->sigmagradPhiFlux[fid] + mesh->sigmaUcrossBFlux[fid];

		}
	}





}

