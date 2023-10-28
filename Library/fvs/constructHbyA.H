#include "mesh.h"

//#include "../Library/fvs/constructHbyA.H"

void constructHbyA(mesh_t *mesh);

void constructHbyA(mesh_t *mesh){

	dfloat *sum;
	dfloat* HbyAf;

	createField(&sum,vector,mesh);
	createField(&HbyAf,surfaceVector,mesh);


	//  construct HbyA and HbyA Fluxes

	for (int cell = 0; cell < mesh->nicells; cell++) {

		for (int i = 0; i < mesh->ndimension; i++) { sum[i]=0.0;}

		for (int face = 0; face < mesh->ncellFaces; face++) {

			int fid 	= face+cell*mesh->ncellFaces;
			int cell_n 	= mesh->faceToNeighborCell[fid];

			for (int dim = 0; dim < mesh->ndimension; dim++) {

				int cell_dim 	=	dim	+ cell*mesh->ndimension;
				int cell_n_dim 	=	dim	+ cell_n*mesh->ndimension;
				int fid_dim 	=	dim	+ fid*mesh->ndimension;

				sum[dim] = sum[dim] + mesh->AAN[fid_dim]*mesh->U[cell_n_dim];
			}

		}

		for (int dim = 0; dim < mesh->ndimension; dim++) {

			int cell_dim = dim+cell*mesh->ndimension;

			mesh->HbyA[cell_dim] = (mesh->BBP[cell_dim] - sum[dim])*mesh->rAP[cell_dim];
		}

	}




	// extrapolate to boundary cells
//		extrapolateToBoundaryCellVector(mesh->HbyA,mesh);




	// constrain HbyA to U for wall boundary condition

	for (dint face = 0; face < mesh->nBoundaryFaces; face++) {

		dint fid   = mesh->boundaryFaces[face];
		dint cell  =  mesh->faceToBoundaryCell[fid];

		for (int dim = 0; dim < mesh->ndimension; dim++) {

			int fid_dim   = dim + fid*mesh->ndimension;
			int cell_dim  = dim + cell*mesh->ndimension;

			mesh->HbyA[cell_dim] = mesh->U[cell_dim];

		}
	}


	// ZeroGradient for Axisymmetric Boundaries

	correctExtrapolationToAxisymmetricBoundaryCellVector(mesh->HbyA,mesh);

	// Find face values

	interpolateCellToFaceLinear(mesh->HbyA,HbyAf,volumeVector,mesh);

	constructFluxofSurfaceVector(HbyAf,mesh->phiHbyA,mesh);

	free(sum);
	free(HbyAf);
}


