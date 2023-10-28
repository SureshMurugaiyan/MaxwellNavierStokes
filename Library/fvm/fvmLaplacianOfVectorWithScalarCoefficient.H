#include "mesh.h"

//#include "../Library/fvm/fvmLaplacianOfVectorWithScalarCoefficient.H"

dfloat fvmLaplacianOfVectorWithScalarCoefficient
(
	dfloat GL,
	dfloat* alpha_face,
	dfloat* phi,
	dint phiName,
	mesh_t *mesh
);

dfloat fvmLaplacianOfVectorWithScalarCoefficient
(
	dfloat GL,
	dfloat* alpha_face,
	dfloat* phi,
	dint phiName,
	mesh_t *mesh
)
{

	/**********************Calculate the Coefficients of Matrix ************************************/

	for (dint face = 0; face < mesh->nInnerAndPeriodicFaces; face++){

		dint fid 	= mesh->innerAndPeriodicFaces[face];
		dint cell 	= mesh->faceToOwnerCell[fid];

		dfloat AP_new  = - GL*mesh->rV[cell]*alpha_face[fid]*mesh->magSf[fid]*mesh->rdelta[fid];
		dfloat AN_new  =  -AP_new;

		for (dint dim = 0; dim < mesh->ndimension; dim++) {

			dint cell_dim = dim + cell*mesh->ndimension;
			dint fid_dim  = dim + fid*mesh->ndimension;

			mesh->AP[cell_dim]	=	mesh->AP[cell_dim] + AP_new;
			mesh->AN[fid_dim]	=	mesh->AN[fid_dim]  + AN_new;
		}

	}


	/**********************Calculate the Coefficients for Boundary Faces *****************************/

	for (dint patch = 0; patch < mesh->npatchs; patch++) {

		dint nPatchFaces = mesh->npatchFaces[patch];

		dfloat* bcValue;
		bcValue = (dfloat*) calloc(nPatchFaces*mesh->ndimension, sizeof(dfloat));

		getBoundaryConditionValue(patch,phiName,bcValue,mesh);

		dint bcType;
		bcType = getBoundaryConditionType(patch,phiName,mesh);

		dint face = mesh->patchStartFaceIndex[patch];

		//==================================================================!
		if (bcType == fixedValue){

			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++) {

				dint fid   = mesh->patchToBoundaryFace[face];
				dint cell  = mesh->faceToOwnerCell[fid];

				for (dint dim = 0; dim < mesh->ndimension; dim++) {

					dint cell_dim        = dim + cell*mesh->ndimension;
					dint PatchFace_dim   = dim + PatchFace*mesh->ndimension;

					dfloat AP_new  = - GL*mesh->rV[cell]*alpha_face[fid]*mesh->magSf[fid]*mesh->rdelta[fid];
					dfloat BP_new  =   AP_new*bcValue[PatchFace_dim];

					mesh->AP[cell_dim]	= mesh->AP[cell_dim] + AP_new;
					mesh->BP[cell_dim]	= mesh->BP[cell_dim] + BP_new;
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

				for (dint dim = 0; dim < mesh->ndimension; dim++) {

					dint cell_dim        = dim + cell*mesh->ndimension;
					dint PatchFace_dim   = dim + PatchFace*mesh->ndimension;

					dfloat coeff  = - GL*mesh->rV[cell]*alpha_face[fid]*mesh->magSf[fid];
					dfloat BP_new =  coeff*bcValue[PatchFace_dim];

					mesh->BP[cell_dim]	= mesh->BP[cell_dim] + BP_new;
				}


				face++;
			}
		}
		//==================================================================!
		if ( bcType == emptyBoundary ){

			// For emptyBoundary no contributions

		}
		//==================================================================!
		free(bcValue);
	}



	return 0.0;
}
