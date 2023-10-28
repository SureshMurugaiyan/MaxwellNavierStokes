#include "mesh.h"

//#include "../Library/fvm/fvmDivergenceOfScalar.H"

dfloat fvmDivergenceOfScalar
(
	dfloat GD,
	dfloat* Phi,
	dint PhiName,
	dfloat* Flux,
	dint interpolationScheme,
	dfloat blendingCoeff,
	mesh_t *mesh
);

dfloat fvmDivergenceOfScalar
(
	dfloat GD,
	dfloat* Phi,
	dint PhiName,
	dfloat* Flux,
	dint interpolationScheme,
	dfloat blendingCoeff,
	mesh_t *mesh
)
{

	/******Calculate the Coefficients of Matrix for innerFaces******/
	dfloat lambda  = 0.0;


	if ((interpolationScheme == linear) || (interpolationScheme == upwind)){

		if (interpolationScheme == linear){
			lambda  = 1.0;
		}
		else if (interpolationScheme == upwind){
			lambda  = 0.0;
		}


		for (dint face = 0; face < mesh->nInnerAndPeriodicFaces; face++){

			dint fid 	= mesh->innerAndPeriodicFaces[face];
			dint cell 	= mesh->faceToOwnerCell[fid];
			dint cell_n	= mesh->faceToNeighborCell[fid];

			dfloat coeff = GD*mesh->rV[cell]*Flux[fid];

			dfloat fxUP = mymax(sgn(Flux[fid]),0.0);
			dfloat fxCD = mesh->fx[fid];


			dfloat AP_new  =  coeff*(1.0-lambda)*fxUP;
			dfloat AN_new  =  coeff*(1.0-lambda)*(1.0-fxUP);
			dfloat BP_new  = -coeff*lambda*( fxCD*Phi[cell] + (1.0-fxCD)*Phi[cell_n]);

			mesh->aP[cell]	=	mesh->aP[cell] + AP_new;
			mesh->bP[cell]	=	mesh->bP[cell] + BP_new;
			mesh->aN[fid]	=	mesh->aN[fid]  + AN_new;


		}

	}
	else{

		if ((blendingCoeff < 0.0) || (blendingCoeff > 1.0)){
			printf ( "\n" );
			printf ( "Fatal error!\n" );
			printf ( " blendingCoeff should be >= 0 and <= 1.0" );
			exit ( 1 );
		}

		if(interpolationScheme == Gamma){
			if ((blendingCoeff < 0.0) || (blendingCoeff > 0.5)){
				printf ( "\n" );
				printf ( "Fatal error!\n" );
				printf ( " For Gamma Scheme blendingCoeff should be >= 0 and <= 0.5" );
				exit ( 1 );
			}
		}




		dfloat* PhiFace; // Face Value of Field
		dfloat *gradPhi; // gradient of field at cell centers

		createField(&PhiFace,surfaceScalar,mesh);
		createField(&gradPhi,volumeVector,mesh);

		interpolateCellToFaceLinear(Phi,PhiFace,volumeScalar, mesh);
		fvdGradientOfScalar(1.0,PhiFace,mesh,gradPhi);

		for (dint face = 0; face < mesh->nInnerAndPeriodicFaces; face++){

			dint fid 	= mesh->innerAndPeriodicFaces[face];
			dint cell 	= mesh->faceToOwnerCell[fid];
			dint cell_n	= mesh->faceToNeighborCell[fid];

			//====================================== Compute Limiter Function
			dfloat phiP =  Phi[cell];
			dfloat phiN =  Phi[cell_n];
			dfloat FaceFlux = Flux[fid];

			dfloat *gradCP;
			dfloat *gradCN;
			dfloat *d;

			gradCP = (dfloat*) calloc(mesh->ndimension, sizeof(dfloat));
			gradCN = (dfloat*) calloc(mesh->ndimension, sizeof(dfloat));
			d	   = (dfloat*) calloc(mesh->ndimension, sizeof(dfloat));


			for (dint dim = 0; dim < mesh->ndimension; dim++) {

				dint cell_dim  = dim + cell*mesh->ndimension;
				dint cell_ndim = dim + cell_n*mesh->ndimension;
				dint fid_dim   = dim + fid*mesh->ndimension;

				gradCP[dim] = gradPhi[cell_dim];
				gradCN[dim] = gradPhi[cell_ndim];
				d[dim] 	    = mesh->d[fid_dim];

			}

			lambda = computeLimiter(FaceFlux,phiP,phiN,gradCP,gradCN,d,interpolationScheme,blendingCoeff,mesh);

			free(gradCP);
			free(gradCN);
			free(d);

			//======================================

			dfloat coeff = GD*mesh->rV[cell]*Flux[fid];

			dfloat fxUP = mymax(sgn(Flux[fid]),0.0);
			dfloat fxCD = mesh->fx[fid];

			dfloat AP_new  =  coeff*(1.0-lambda)*fxUP;
			dfloat AN_new  =  coeff*(1.0-lambda)*(1.0-fxUP);
			dfloat BP_new  = -coeff*lambda*( fxCD*Phi[cell] + (1.0-fxCD)*Phi[cell_n]);

			mesh->aP[cell]	=	mesh->aP[cell] + AP_new;
			mesh->bP[cell]	=	mesh->bP[cell] + BP_new;
			mesh->aN[fid]	=	mesh->aN[fid]  + AN_new;

			//======================================

		}

		free(PhiFace);
		free(gradPhi);

	}




	/**********************Calculate the Coefficients for Boundary Faces *****************************/

	for (dint patch = 0; patch < mesh->npatchs; patch++) {

		dint nPatchFaces = mesh->npatchFaces[patch];


		dfloat* bcValue;
		bcValue = (dfloat*) calloc(nPatchFaces, sizeof(dfloat));

		getBoundaryConditionValue(patch,PhiName,bcValue, mesh);

		dint bcType;
		bcType = getBoundaryConditionType(patch,PhiName,mesh);

		dint face = mesh->patchStartFaceIndex[patch];

		//==================================================================!
		if (bcType == fixedValue){

			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++) {

				dint fid  = mesh->patchToBoundaryFace[face];
				dint cell = mesh->faceToOwnerCell[fid];

				dfloat coeff   =  -GD*mesh->rV[cell]*Flux[fid];

				dfloat R_new   =    coeff*bcValue[PatchFace];

				mesh->bP[cell]	=	mesh->bP[cell]  + R_new;

				face++;

			}
		}

		//==================================================================!
		if (bcType == zeroNormalGradient){

			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++) {

				dint fid  = mesh->patchToBoundaryFace[face];
				dint cell = mesh->faceToOwnerCell[fid];

				dfloat AP_new   =  GD*mesh->rV[cell]*Flux[fid];

				mesh->aP[cell]	=	mesh->aP[cell]  + AP_new;

				face++;

			}

		}
		//==================================================================!

		if (bcType == fixedNormalGradient){

			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++) {

				dint fid  = mesh->patchToBoundaryFace[face];
				dint cell = mesh->faceToOwnerCell[fid];

				dfloat coeff   =  GD*mesh->rV[cell]*Flux[fid];

				mesh->aP[cell]	=	mesh->aP[cell]  + coeff;
				mesh->bP[cell]	=	mesh->bP[cell]  - coeff*bcValue[PatchFace]*mesh->magdf[fid];


				face++;

			}

		}
		//==================================================================!

		free(bcValue);

	}



	return 0.0;

}
