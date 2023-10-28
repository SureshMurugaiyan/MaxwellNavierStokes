#include "mesh.h"

//#include "../Library/fvs/correctExtrapolationToAxisymmetricBoundaryCellScalar.H"


void correctExtrapolationToAxisymmetricBoundaryCellScalar(dfloat* phiCell,mesh_t *mesh);

void correctExtrapolationToAxisymmetricBoundaryCellScalar(dfloat* phiCell,mesh_t *mesh){


	// Zero Gradient for Axisymmetric boundary faces

	for (dint face = 0; face < mesh->nBoundaryFaces; face++) {

		dint fid     = mesh->boundaryFaces[face];
		dint Ocell   = mesh->faceToOwnerCell[fid];
		dint bcell   = mesh->faceToBoundaryCell[fid];

		dint lface = fid % mesh->ncellFaces; // local face

		if(lface==mesh->leftFace)   {

			phiCell[bcell] = phiCell[Ocell];

		}

		if(lface==mesh->rightFace)   {

			phiCell[bcell] = phiCell[Ocell];
		}



	}


//	strcpy(mesh->variableName, "phiCell");
//	printVolumeScalar(phiCell, mesh);
}
