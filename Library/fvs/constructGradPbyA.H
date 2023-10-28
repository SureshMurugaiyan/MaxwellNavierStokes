#include "mesh.h"

//#include "../Library/fvs/constructGradPbyA.H"

void constructGradPbyA(mesh_t *mesh);

void constructGradPbyA(mesh_t *mesh){

	dfloat *pf;    // cell face values for dynamic pressure
	dfloat *gradP; // cell center Pressure gradient
	dfloat* Sfm;
	dfloat* magSfm;

	createField(&pf,surfaceScalar,mesh);
	createField(&gradP,volumeVector,mesh);
	createField(&Sfm,surfaceVector,mesh);
	createField(&magSfm,surfaceScalar,mesh);

	interpolateCellToFaceLinear(mesh->p,pf,volumeScalar,mesh);

	fvdGradientOfScalar(1.0, pf, mesh, gradP);

	for (int cell = 0; cell < mesh->nicells; cell++) {
		for (int dim = 0; dim < mesh->ndimension; dim++) {
			int cell_dim = dim+ cell*mesh->ndimension;

			mesh->gradPbyA[cell_dim] = gradP[cell_dim]*mesh->rAP[cell_dim];
		}
	}



	// Flux of gradPmbyA

	computeModifiedFaceArea(mesh->rAPf, Sfm,magSfm,mesh);


	for (int cell = 0; cell < mesh->nicells; cell++) {

		for (int face = 0; face < mesh->ncellFaces; face++) {

			dint fid 	= face + cell*mesh->ncellFaces;
			int cell_n	=  mesh->faceToNeighborCell[fid];

			mesh->phigradPbyA[fid] = (mesh->p[cell_n]-mesh->p[cell])*magSfm[fid]*mesh->rdelta[fid];
		}

	}


	free(pf);
	free(gradP);
	free(Sfm);
	free(magSfm);
}


