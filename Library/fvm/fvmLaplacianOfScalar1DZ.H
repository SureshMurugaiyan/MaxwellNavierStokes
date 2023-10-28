#include "mesh.h"

//#include "../Library/fvm/fvmLaplacianOfScalar1DZ.H"

dfloat fvmLaplacianOfScalar1DZ
(
	dfloat GL,
	dfloat* phi,
	dint phiName,
	mesh_t *mesh
);

dfloat fvmLaplacianOfScalar1DZ
(
	dfloat GL,
	dfloat* phi,
	dint phiName,
	mesh_t *mesh
)
{


	/**********************Calculate the Coefficients of Matrix ************************************/

	for (dint face = 0; face < mesh->nInnerAndPeriodicFaces; face++){

		dint fid 	= mesh->innerAndPeriodicFaces[face];
		dint cell 	= mesh->faceToOwnerCell[fid];

		int lface = fid % mesh->ncellFaces;

		if((lface == mesh->topFace) ||  (lface==mesh->bottomFace)){

			dfloat aP_new  = - GL*mesh->rV[cell]*mesh->magSf[fid]*mesh->rmagd[fid];
			dfloat aN_new  = - aP_new;

			mesh->aP[cell]	=	mesh->aP[cell] + aP_new;
			mesh->aN[fid]	=	mesh->aN[fid]  + aN_new;


		}


	}


	/**********************Calculate the Coefficients for Boundary Faces *****************************/

	for (dint patch = 0; patch < mesh->npatchs; patch++) {

		dint nPatchFaces = mesh->npatchFaces[patch];

		dfloat* bcValue;
		bcValue = (dfloat*) calloc(nPatchFaces, sizeof(dfloat));

		getBoundaryConditionValue(patch,phiName,bcValue,mesh);

		dint bcType;
		bcType = getBoundaryConditionType(patch,phiName,mesh);

		dint face = mesh->patchStartFaceIndex[patch];

		//==================================================================!
		if (bcType == fixedValue){

			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++) {

				dint fid   = mesh->patchToBoundaryFace[face];
				dint cell  = mesh->faceToOwnerCell[fid];

				int lface = fid % mesh->ncellFaces;

				if((lface == mesh->topFace) ||  (lface==mesh->bottomFace)){

					dfloat AP_new  = - GL*mesh->rV[cell]*mesh->magSf[fid]*mesh->rmagdf[fid];

					dfloat R_new   =    AP_new*bcValue[PatchFace];

					mesh->aP[cell]	=	mesh->aP[cell] + AP_new;
					mesh->bP[cell]	=	mesh->bP[cell]  + R_new;


				}


				face++;
			}

		}
		//==================================================================!
		if (bcType == zeroNormalGradient){

			// For zero gradient term all terms vanish

		}
		//==================================================================!

		if (bcType == fixedNormalGradient){

			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++) {

				dint fid   = mesh->patchToBoundaryFace[face];
				dint cell  = mesh->faceToOwnerCell[fid];

				int lface = fid % mesh->ncellFaces;

				if((lface == mesh->topFace) ||  (lface==mesh->bottomFace)){

					dfloat coeff   = - GL*mesh->rV[cell]*mesh->magSf[fid]*bcValue[PatchFace];

					mesh->bP[cell]	= mesh->bP[cell] + coeff;


				}

				face++;
			}

		}
		//==================================================================!
		free(bcValue);
	}



	return 0.0;
}
