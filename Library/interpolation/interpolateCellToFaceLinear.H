#include "mesh.h"

//#include "../Library/interpolation/interpolateCellToFaceLinear.H"

void interpolateCellToFaceLinear
(
	dfloat* phiCell,
	dfloat* phiFace,
	dint fieldType,
	mesh_t*	mesh
);

void interpolateCellToFaceLinear
(
	dfloat* phiCell,
	dfloat* phiFace,
	dint fieldType,
	mesh_t*	mesh
)
{

	dint ndimension   = getNdimension(fieldType);

	for (int face = 0; face < mesh->nInnerAndPeriodicFaces; face++) {

		int fid 	=  mesh->innerAndPeriodicFaces[face];
		int cell	=  mesh->faceToOwnerCell[fid];
		int cell_n	=  mesh->faceToNeighborCell[fid];

		for (int dim = 0; dim < ndimension; dim++) {

			int fid_dim   = dim + fid*ndimension;
			int cell_dim  = dim + cell*ndimension;
			int cell_ndim = dim + cell_n*ndimension;

			phiFace[fid_dim] = mesh->fx[fid]*phiCell[cell_dim] +
							   (1-mesh->fx[fid])*phiCell[cell_ndim];
		}
	}



	for (int patch = 0; patch < mesh->npatchs; patch++) {

		int nPatchFaces = mesh->npatchFaces[patch];

		int npatch  = mesh->periodicNeighborPatch[patch];
		int face    = mesh->patchStartFaceIndex[patch];

		if(npatch == -1){

			for (int PatchFace = 0; PatchFace < nPatchFaces; PatchFace++) {

				int fid   		 = mesh->patchToBoundaryFace[face];
				int Boundarycell = mesh->faceToBoundaryCell[fid];

				for (int dim = 0; dim < ndimension; dim++) {

					int Boundarycell_dim  = dim + Boundarycell*ndimension;
					int fid_dim           = dim + fid*ndimension;

					phiFace[fid_dim] = phiCell[Boundarycell_dim];
				}

				face++;
			}

		}
	}


}

