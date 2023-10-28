#include "mesh.h"

//#include "../Library/fvm/fvmLaplacianOfScalarWithScalarCoefficient.H"

dfloat fvmLaplacianOfScalarWithScalarCoefficient
(
	dfloat GL,
	dfloat* alpha_face,
	dfloat* phi,
	dint phiName,
	mesh_t *mesh
);

dfloat fvmLaplacianOfScalarWithScalarCoefficient
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

		dfloat aP_new  = - GL*mesh->rV[cell]*alpha_face[fid]*mesh->magSf[fid]*mesh->rdelta[fid];
		dfloat aN_new  =  -aP_new;

		mesh->aP[cell]	=	mesh->aP[cell] + aP_new;
		mesh->aN[fid]	=	mesh->aN[fid]  + aN_new;

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

				dfloat aP_new  = - GL*mesh->rV[cell]*alpha_face[fid]*mesh->magSf[fid]*mesh->rdelta[fid];
				dfloat r_new   =   aP_new*bcValue[PatchFace];

				mesh->aP[cell]	=	mesh->aP[cell] + aP_new;
				mesh->bP[cell]	=	mesh->bP[cell] + r_new;

				face++;
			}

		}
		//==================================================================!
		if ( bcType == zeroNormalGradient ){

			// For zero gradient term all terms vanish

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
