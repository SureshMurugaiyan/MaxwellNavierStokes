#include "mesh.h"

//#include "../Library/fvs/updateBoundaryValues.H"


void updateBoundaryValues(dfloat* phi, int phiName, mesh_t *mesh);

void updateBoundaryValues(dfloat* phi, int phiName, mesh_t *mesh){


	dint fieldType    = getFieldType(phiName);
	dint ndimension   = getNdimension(fieldType);

	for (int patch = 0; patch < mesh->npatchs; patch++) {

		int nPatchFaces = mesh->npatchFaces[patch];

		dfloat* bcValue;
		bcValue = (dfloat*) calloc(nPatchFaces*ndimension, sizeof(dfloat));

		getBoundaryConditionValue(patch,phiName,bcValue, mesh);

		int bcType;
		bcType = getBoundaryConditionType(patch,phiName,mesh);

		int face = mesh->patchStartFaceIndex[patch];

		if (bcType == fixedValue){

			for (int PatchFace = 0; PatchFace < nPatchFaces; PatchFace++) {

				int fid   		 = mesh->patchToBoundaryFace[face];
				int Boundarycell = mesh->faceToBoundaryCell[fid];

				for (int dim = 0; dim < ndimension; dim++) {

					int Boundarycell_dim  = dim + Boundarycell*ndimension;
					int PatchFace_dim     = dim + PatchFace*ndimension;

					phi[Boundarycell_dim] = bcValue[PatchFace_dim];

				}
				face++;
			}
		}
		//==================================================================!
		if (bcType == zeroNormalGradient){

			for (int PatchFace = 0; PatchFace < nPatchFaces; PatchFace++) {

				int fid   		 = mesh->patchToBoundaryFace[face];
				int Boundarycell = mesh->faceToBoundaryCell[fid];
				int Ownercell 	 = mesh->faceToOwnerCell[fid];

				for (int dim = 0; dim < ndimension; dim++) {

					int Bdarycell_dim  = dim + Boundarycell*ndimension;
					int Ownercell_dim  = dim + Ownercell*ndimension;

					phi[Bdarycell_dim] = phi[Ownercell_dim];

				}

				face++;
			}
		}
		//==================================================================!
		if (bcType == periodic){

			// No need to update boundary values, all are interior values
			// These boundary values are not used in calculations
			// For now copying the values using zeroNormalGradient


			for (int PatchFace = 0; PatchFace < nPatchFaces; PatchFace++) {

				int fid   		 = mesh->patchToBoundaryFace[face];
				int Boundarycell = mesh->faceToBoundaryCell[fid];
				int Ownercell 	 = mesh->faceToOwnerCell[fid];

				for (int dim = 0; dim < ndimension; dim++) {

					int Bdarycell_dim  = dim + Boundarycell*ndimension;
					int Ownercell_dim  = dim + Ownercell*ndimension;

					phi[Bdarycell_dim] = phi[Ownercell_dim];

				}

				face++;
			}


		}
		//==================================================================!

		if (bcType == fixedNormalGradient){

			for (int PatchFace = 0; PatchFace < nPatchFaces; PatchFace++) {

				int fid   		 = mesh->patchToBoundaryFace[face];
				int Boundarycell = mesh->faceToBoundaryCell[fid];
				int Ownercell 	 = mesh->faceToOwnerCell[fid];

				for (int dim = 0; dim < ndimension; dim++) {

					int Boundarycell_dim  = dim + Boundarycell*ndimension;
					int PatchFace_dim     = dim + PatchFace*ndimension;
					int Ownercell_dim     = dim + Ownercell*ndimension;


					phi[Boundarycell_dim] = phi[Ownercell_dim] + bcValue[PatchFace_dim]*mesh->magdf[fid];

				}
				face++;
			}
		}

		//==================================================================!

		if (bcType == emptyBoundary){

			// Don't update the boundary cells for empty boundary
		}



		//==================================================================!
		free(bcValue);

	}


}
