#include "mesh.h"

//#include "../Library/fvs/correctExtrapolationToAxisymmetricBoundaryCellVector.H"

void correctExtrapolationToAxisymmetricBoundaryCellVector(dfloat* phiCell,mesh_t *mesh);

void correctExtrapolationToAxisymmetricBoundaryCellVector(dfloat* phiCell,mesh_t *mesh){

	// Zero Gradient for Axisymmetric boundary faces

	for (dint face = 0; face < mesh->nBoundaryFaces; face++) {

		dint fid     = mesh->boundaryFaces[face];
		dint Ocell   = mesh->faceToOwnerCell[fid];
		dint Bcell   = mesh->faceToBoundaryCell[fid];

		dint lface = fid % mesh->ncellFaces; // local face

		if(lface==mesh->leftFace)   {

			for (int dim = 0; dim < mesh->ndimension; dim++) {

				int Ocell_dim = dim + Ocell*mesh->ndimension;
				int Bcell_dim = dim + Bcell*mesh->ndimension;

				phiCell[Bcell_dim] = phiCell[Ocell_dim];

			}

		}

		if(lface==mesh->rightFace)   {

			for (int dim = 0; dim < mesh->ndimension; dim++) {

				int Ocell_dim = dim + Ocell*mesh->ndimension;
				int Bcell_dim = dim + Bcell*mesh->ndimension;

				phiCell[Bcell_dim] = phiCell[Ocell_dim];

			}

		}



	}


//	strcpy(mesh->variableName, "phiCell");
//	printVolumeVector(phiCell, mesh);

}
