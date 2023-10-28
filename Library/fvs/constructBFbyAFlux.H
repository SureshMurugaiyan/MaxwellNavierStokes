#include "mesh.h"

//#include "../Library/fvs/constructBFbyAFlux.H"

void constructBFbyAFlux(mesh_t *mesh);

void constructBFbyAFlux(mesh_t *mesh){


	for (int cell = 0; cell < mesh->nicells; cell++) {
		for (int dim = 0; dim < mesh->ndimension; dim++) {
			int cell_dim = dim+ cell*mesh->ndimension;

			mesh->BF1byA[cell_dim] = mesh->BF1[cell_dim]*mesh->rAP[cell_dim];
		}
	}


	if(mesh->useOriginalGravityFormulation	== true){


		dfloat* BF1UbyAf;
		createField(&BF1UbyAf,surfaceVector,mesh);


		// Flux of BF1byA
		interpolateCellToFaceLinear(mesh->BF1byA,BF1UbyAf,volumeVector,mesh);
		constructFluxofSurfaceVector(BF1UbyAf,mesh->phiBF1byA,mesh);


		free(BF1UbyAf);

	}


	if(mesh->useModifiedGravityFormulation	== true){


		dfloat* Sfm;
		dfloat* magSfm;

		createField(&Sfm,surfaceVector,mesh);
		createField(&magSfm,surfaceScalar,mesh);



		// Flux of BF1byA

		computeModifiedFaceArea(mesh->rAPf, Sfm,magSfm,mesh);

		dfloat oneOverFrsqr = 1.0/(mesh->Fr*mesh->Fr);

		for (int cell = 0; cell < mesh->nicells; cell++) {

			for (int face = 0; face < mesh->ncellFaces; face++) {

				dint fid 	= face + cell*mesh->ncellFaces;
				int cell_n	=  mesh->faceToNeighborCell[fid];

				mesh->phiBF1byA[fid] = oneOverFrsqr*mesh->gDotXf[fid]*(mesh->rho[cell_n]-mesh->rho[cell])
						                                             *mesh->rdelta[fid]*magSfm[fid];
			}

		}


		free(Sfm);
		free(magSfm);


	}



}
