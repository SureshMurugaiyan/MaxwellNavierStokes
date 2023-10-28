#include "mesh.h"

//#include "../Library/fvm/fvmLaplacianOfScalarWithVectorCoefficient.H"

dfloat fvmLaplacianOfScalarWithVectorCoefficient
(
	dfloat GL,
	dfloat* alpha_f,
	dfloat* phi,
	dint phiName,
	mesh_t *mesh
);

dfloat fvmLaplacianOfScalarWithVectorCoefficient
(
	dfloat GL,
	dfloat* alpha_f,
	dfloat* phi,
	dint phiName,
	mesh_t *mesh
)
{

	// Compute the modified face area
	dfloat* Sfm;
	dfloat* magSfm;

	createField(&Sfm,surfaceVector,mesh);
	createField(&magSfm,surfaceScalar,mesh);

	computeModifiedFaceArea(alpha_f, Sfm,magSfm,mesh);
	/**********************Calculate the Coefficients of Matrix ************************************/

	for (dint face = 0; face < mesh->nInnerAndPeriodicFaces; face++){

		dint fid 	= mesh->innerAndPeriodicFaces[face];
		dint cell 	= mesh->faceToOwnerCell[fid];


		dfloat aP_new  = - GL*mesh->rV[cell]*magSfm[fid]*mesh->rmagd[fid];
		dfloat aN_new  = - aP_new;

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

				dfloat AP_new  = - GL*mesh->rV[cell]*magSfm[fid]*mesh->rmagdf[fid];

				dfloat R_new   =    AP_new*bcValue[PatchFace];

				mesh->aP[cell]	=	mesh->aP[cell] + AP_new;
				mesh->bP[cell]	=	mesh->bP[cell]  + R_new;

				face++;
			}

		}
		//==================================================================!
		if (bcType == zeroNormalGradient){

			// For zero gradient term all terms vanish

		}
		//==================================================================!
		if (bcType == totalPressureInletPBC){

			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++) {

				dint fid   = mesh->patchToBoundaryFace[face];
				dint cell  = mesh->faceToOwnerCell[fid];


				dint idx	= mesh->dim_x + cell*mesh->ndimension;
				dint idy	= mesh->dim_y + cell*mesh->ndimension;
				dint idz	= mesh->dim_z + cell*mesh->ndimension;

				dfloat magUsqr = 0.0;

				dfloat a = mesh->U[idx];
				dfloat b = mesh->U[idy];
				dfloat c = mesh->U[idz];

				magUsqr = (a*a +b*b +c*c);

				dfloat pValue = 1.0 - 0.5*magUsqr;

				dfloat AP_new  = - GL*mesh->rV[cell]*magSfm[fid]*mesh->rmagdf[fid];

				dfloat R_new   =    AP_new*pValue;

				mesh->aP[cell]	=	mesh->aP[cell] + AP_new;
				mesh->bP[cell]	=	mesh->bP[cell]  + R_new;

				face++;
			}

		}
		//==================================================================!
		if (bcType == totalPressureOutletPBC){

			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++) {

				dint fid   = mesh->patchToBoundaryFace[face];
				dint cell  = mesh->faceToOwnerCell[fid];


				dint idx	= mesh->dim_x + cell*mesh->ndimension;
				dint idy	= mesh->dim_y + cell*mesh->ndimension;
				dint idz	= mesh->dim_z + cell*mesh->ndimension;

				dfloat magUsqr = 0.0;

				dfloat a = mesh->U[idx];
				dfloat b = mesh->U[idy];
				dfloat c = mesh->U[idz];

				magUsqr = (a*a +b*b +c*c);

				dfloat pValue = 1.0 ;

				dfloat AP_new  = - GL*mesh->rV[cell]*magSfm[fid]*mesh->rmagdf[fid];

				dfloat R_new   =    AP_new*pValue;

				mesh->aP[cell]	=	mesh->aP[cell] + AP_new;
				mesh->bP[cell]	=	mesh->bP[cell]  + R_new;

				face++;
			}

		}
		//==================================================================!

		if (bcType == fixedNormalGradient){

			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++) {

				dint fid   = mesh->patchToBoundaryFace[face];
				dint cell  = mesh->faceToOwnerCell[fid];

				dfloat coeff   = - GL*mesh->rV[cell]*magSfm[fid]*bcValue[PatchFace];

				mesh->bP[cell]	= mesh->bP[cell] + coeff;

				face++;
			}

		}
		//==================================================================!
		free(bcValue);
	}

	free(Sfm);
	free(magSfm);

	return 0.0;
}
