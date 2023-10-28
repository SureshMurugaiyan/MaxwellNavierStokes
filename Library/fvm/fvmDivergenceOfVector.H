#include "mesh.h"

//#include "../Library/fvm/fvmDivergenceOfVector.H"

dfloat fvmDivergenceOfVector
(
	dfloat GD,
	dfloat* Phi,
	dint PhiName,
	dfloat* Flux,
	dint interpolationScheme,
	dfloat blendingCoeff,
	mesh_t *mesh
);

dfloat fvmDivergenceOfVector
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

			for (dint dim = 0; dim < mesh->ndimension; dim++) {

				dint cell_dim   = dim + cell*mesh->ndimension;
				dint cell_n_dim = dim + cell_n*mesh->ndimension;
				dint fid_dim    = dim + fid*mesh->ndimension;


				dfloat coeff = GD*mesh->rV[cell]*Flux[fid];

				dfloat fxUP = mymax(sgn(Flux[fid]),0.0);
				dfloat fxCD = mesh->fx[fid];


				dfloat AP_new  =  coeff*(1.0-lambda)*fxUP;
				dfloat AN_new  =  coeff*(1.0-lambda)*(1.0-fxUP);
				dfloat BP_new  = -coeff*lambda*( fxCD*Phi[cell_dim] + (1.0-fxCD)*Phi[cell_n_dim]);

				mesh->AP[cell_dim]	=	mesh->AP[cell_dim] + AP_new;
				mesh->BP[cell_dim]	=	mesh->BP[cell_dim] + BP_new;
				mesh->AN[fid_dim]	=	mesh->AN[fid_dim]  + AN_new;

			}




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


		dfloat *Phix; // Extract x-compoent of the field
		dfloat *Phiy; // Extract y-compoent of the field
		dfloat *Phiz; // Extract z-compoent of the field

		createField(&Phix,volumeScalar,mesh);
		createField(&Phiy,volumeScalar,mesh);
		createField(&Phiz,volumeScalar,mesh);


		for(int cell=0;cell<mesh->ncells;++cell){

			dint dimx = 0; // dim = 0 gives x component
			dint dimy = 1; // dim = 1 gives y component
			dint dimz = 2; // dim = 2 gives z component

			dint idx	= dimx+cell*mesh->ndimension;
			dint idy	= dimy+cell*mesh->ndimension;
			dint idz	= dimz+cell*mesh->ndimension;

			Phix[cell]	= Phi[idx];
			Phiy[cell]	= Phi[idy];
			Phiz[cell]	= Phi[idz];


		}


		dfloat* PhiFacex; // Face Value of Field
		dfloat* PhiFacey; // Face Value of Field
		dfloat* PhiFacez; // Face Value of Field

		dfloat *gradPhix; // gradient of field at cell centers
		dfloat *gradPhiy; // gradient of field at cell centers
		dfloat *gradPhiz; // gradient of field at cell centers

		createField(&PhiFacex,surfaceScalar,mesh);
		createField(&PhiFacey,surfaceScalar,mesh);
		createField(&PhiFacez,surfaceScalar,mesh);

		createField(&gradPhix,volumeVector,mesh);
		createField(&gradPhiy,volumeVector,mesh);
		createField(&gradPhiz,volumeVector,mesh);

		interpolateCellToFaceLinear(Phix,PhiFacex,volumeScalar, mesh);
		interpolateCellToFaceLinear(Phiy,PhiFacey,volumeScalar, mesh);
		interpolateCellToFaceLinear(Phiz,PhiFacez,volumeScalar, mesh);

		fvdGradientOfScalar(1.0,PhiFacex,mesh,gradPhix);
		fvdGradientOfScalar(1.0,PhiFacey,mesh,gradPhiy);
		fvdGradientOfScalar(1.0,PhiFacez,mesh,gradPhiz);


		for (dint face = 0; face < mesh->nInnerAndPeriodicFaces; face++){

			dint fid 	= mesh->innerAndPeriodicFaces[face];
			dint cell 	= mesh->faceToOwnerCell[fid];
			dint cell_n	= mesh->faceToNeighborCell[fid];

			dint dimx = 0; // dim = 0 gives x component
			dint dimy = 1; // dim = 1 gives y component
			dint dimz = 2; // dim = 2 gives z component

			dint idx	= dimx+cell*mesh->ndimension;
			dint idy	= dimy+cell*mesh->ndimension;
			dint idz	= dimz+cell*mesh->ndimension;

			dint fidx	= dimx+fid*mesh->ndimension;
			dint fidy	= dimy+fid*mesh->ndimension;
			dint fidz	= dimz+fid*mesh->ndimension;

			//====================================== Compute Limiter Function
			dfloat phiPx =  Phix[cell];
			dfloat phiPy =  Phiy[cell];
			dfloat phiPz =  Phiz[cell];

			dfloat phiNx =  Phix[cell_n];
			dfloat phiNy =  Phiy[cell_n];
			dfloat phiNz =  Phiz[cell_n];

			dfloat FaceFlux = Flux[fid];

			dfloat *gradCPx;
			dfloat *gradCPy;
			dfloat *gradCPz;

			dfloat *gradCNx;
			dfloat *gradCNy;
			dfloat *gradCNz;

			dfloat *d;

			gradCPx = (dfloat*) calloc(mesh->ndimension, sizeof(dfloat));
			gradCPy = (dfloat*) calloc(mesh->ndimension, sizeof(dfloat));
			gradCPz = (dfloat*) calloc(mesh->ndimension, sizeof(dfloat));

			gradCNx = (dfloat*) calloc(mesh->ndimension, sizeof(dfloat));
			gradCNy = (dfloat*) calloc(mesh->ndimension, sizeof(dfloat));
			gradCNz = (dfloat*) calloc(mesh->ndimension, sizeof(dfloat));

			d	   = (dfloat*) calloc(mesh->ndimension, sizeof(dfloat));


			for (dint dim = 0; dim < mesh->ndimension; dim++) {

				dint cell_dim  = dim + cell*mesh->ndimension;
				dint cell_ndim = dim + cell_n*mesh->ndimension;
				dint fid_dim   = dim + fid*mesh->ndimension;

				gradCPx[dim] = gradPhix[cell_dim];
				gradCPy[dim] = gradPhiy[cell_dim];
				gradCPz[dim] = gradPhiz[cell_dim];


				gradCNx[dim] = gradPhix[cell_ndim];
				gradCNy[dim] = gradPhiy[cell_ndim];
				gradCNz[dim] = gradPhiz[cell_ndim];

				d[dim] 	    = mesh->d[fid_dim];

			}

			dfloat lambdax = 0.0;
			dfloat lambday = 0.0;
			dfloat lambdaz = 0.0;


			lambdax = computeLimiter(FaceFlux,phiPx,phiNx,gradCPx,gradCNx,d,interpolationScheme,blendingCoeff,mesh);
			lambday = computeLimiter(FaceFlux,phiPy,phiNy,gradCPy,gradCNy,d,interpolationScheme,blendingCoeff,mesh);
			lambdaz = computeLimiter(FaceFlux,phiPz,phiNz,gradCPz,gradCNz,d,interpolationScheme,blendingCoeff,mesh);


			free(gradCPx);
			free(gradCPy);
			free(gradCPz);

			free(gradCNx);
			free(gradCNy);
			free(gradCNz);

			free(d);

			//======================================

			dfloat coeff = GD*mesh->rV[cell]*Flux[fid];

			dfloat fxUP = mymax(sgn(Flux[fid]),0.0);
			dfloat fxCD = mesh->fx[fid];

			dfloat APx_new  =  coeff*(1.0-lambdax)*fxUP;
			dfloat APy_new  =  coeff*(1.0-lambday)*fxUP;
			dfloat APz_new  =  coeff*(1.0-lambdaz)*fxUP;

			dfloat ANx_new  =  coeff*(1.0-lambdax)*(1.0-fxUP);
			dfloat ANy_new  =  coeff*(1.0-lambday)*(1.0-fxUP);
			dfloat ANz_new  =  coeff*(1.0-lambdaz)*(1.0-fxUP);

			dfloat BPx_new  = -coeff*lambdax*( fxCD*Phix[cell] + (1.0-fxCD)*Phix[cell_n]);
			dfloat BPy_new  = -coeff*lambday*( fxCD*Phiy[cell] + (1.0-fxCD)*Phiy[cell_n]);
			dfloat BPz_new  = -coeff*lambdaz*( fxCD*Phiz[cell] + (1.0-fxCD)*Phiz[cell_n]);

			mesh->AP[idx]	=	mesh->AP[idx] + APx_new;
			mesh->AP[idy]	=	mesh->AP[idy] + APy_new;
			mesh->AP[idz]	=	mesh->AP[idz] + APz_new;


			mesh->BP[idx]	=	mesh->BP[idx] + BPx_new;
			mesh->BP[idy]	=	mesh->BP[idy] + BPy_new;
			mesh->BP[idz]	=	mesh->BP[idz] + BPz_new;


			mesh->AN[fidx]	=	mesh->AN[fidx]  + ANx_new;
			mesh->AN[fidy]	=	mesh->AN[fidy]  + ANy_new;
			mesh->AN[fidz]	=	mesh->AN[fidz]  + ANz_new;


			//======================================

		}

		free(PhiFacex);
		free(PhiFacey);
		free(PhiFacez);

		free(gradPhix);
		free(gradPhiy);
		free(gradPhiz);


		free(Phix);
		free(Phiy);
		free(Phiz);


	}




	/**********************Calculate the Coefficients for Boundary Faces *****************************/

	for (dint patch = 0; patch < mesh->npatchs; patch++) {

		dint nPatchFaces = mesh->npatchFaces[patch];


		dfloat* bcValue;
		bcValue = (dfloat*) calloc(nPatchFaces*mesh->ndimension, sizeof(dfloat));

		getBoundaryConditionValue(patch,PhiName,bcValue, mesh);

		dint bcType;
		bcType = getBoundaryConditionType(patch,PhiName,mesh);

		dint face = mesh->patchStartFaceIndex[patch];

		//==================================================================!
		if (bcType == fixedValue){
			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++) {

				dint fid  = mesh->patchToBoundaryFace[face];
				dint cell = mesh->faceToOwnerCell[fid];

				for (dint dim = 0; dim < mesh->ndimension; dim++) {

					dint cell_dim      = dim + cell*mesh->ndimension;
					dint PatchFace_dim = dim + PatchFace*mesh->ndimension;

					dfloat coeff   =  -GD*mesh->rV[cell]*Flux[fid];

					dfloat R_new   =    coeff*bcValue[PatchFace_dim];

					mesh->BP[cell_dim]	=	mesh->BP[cell_dim]  + R_new;
				}

				face++;

			}
		}

		//==================================================================!
		if (bcType == zeroNormalGradient ){

			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++) {

				dint fid  = mesh->patchToBoundaryFace[face];
				dint cell = mesh->faceToOwnerCell[fid];

				dfloat AP_new   =  GD*mesh->rV[cell]*Flux[fid];

				for (dint dim = 0; dim < mesh->ndimension; dim++) {

						dint cell_dim  = dim+cell*mesh->ndimension;

						mesh->AP[cell_dim]	=	mesh->AP[cell_dim]  + AP_new;
				}



				face++;

			}

		}
	//==================================================================!
		if (bcType == fixedNormalGradient ){

			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++) {

				dint fid  = mesh->patchToBoundaryFace[face];
				dint cell = mesh->faceToOwnerCell[fid];

				dfloat coeff   =  GD*mesh->rV[cell]*Flux[fid];

				for (dint dim = 0; dim < mesh->ndimension; dim++) {

						dint cell_dim  = dim+cell*mesh->ndimension;
						dint PatchFace_dim = dim + PatchFace*mesh->ndimension;

						mesh->AP[cell_dim]	=	mesh->AP[cell_dim]  + coeff;
						mesh->BP[cell_dim]	=	mesh->BP[cell_dim]  - coeff*bcValue[PatchFace_dim]*mesh->magdf[fid];
				}



				face++;

			}
		}
	//==================================================================!
		free(bcValue);

	}






	return 0.0;
}

