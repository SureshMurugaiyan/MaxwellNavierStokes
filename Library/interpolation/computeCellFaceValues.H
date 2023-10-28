#include "mesh.h"

//#include "../Library/interpolation/computeCellFaceValues.H"


void computeCellFaceValues
(
	dfloat* phiCell,
	dfloat* phiFace,
	dint 	fieldType,
	dfloat* Flux,
	dint 	interpolationScheme,
	dfloat 	blendingCoeff,
	mesh_t*	mesh
);

void computeCellFaceValues
(
	dfloat* phiCell,
	dfloat* phiFace,
	dint 	fieldType,
	dfloat* Flux,
	dint 	interpolationScheme,
	dfloat 	blendingCoeff,
	mesh_t*	mesh
)
{

	dint ndimension   = getNdimension(fieldType);

	if(interpolationScheme == linear){
		// Central Scheme
		for (dint face = 0; face < mesh->nInnerAndPeriodicFaces; face++){

			dint fid 	= mesh->innerAndPeriodicFaces[face];
			dint cell 	= mesh->faceToOwnerCell[fid];
			dint cell_n	= mesh->faceToNeighborCell[fid];

			dfloat fxUP = mymax(sgn(Flux[fid]),0.0);
			dfloat fxCD = mesh->fx[fid];

			dfloat lambda  = 1.0;

			for (int dim = 0; dim < ndimension; dim++) {

				dint fid_dim   = dim + fid*ndimension;
				dint cell_dim  = dim + cell*ndimension;
				dint cell_ndim = dim + cell_n*ndimension;

				dfloat phiP =  phiCell[cell_dim];
				dfloat phiN =  phiCell[cell_ndim];

				dfloat UP =  fxUP*phiP + (1-fxUP)*phiN;
				dfloat CD =  fxCD*phiP + (1-fxCD)*phiN;

				phiFace[fid_dim] = (1.0-lambda)*UP + lambda*CD;

			}
		}
	}



	if(interpolationScheme == upwind){
		// UPWIND SCHEME
		for (dint face = 0; face < mesh->nInnerAndPeriodicFaces; face++){

			dint fid 	= mesh->innerAndPeriodicFaces[face];
			dint cell 	= mesh->faceToOwnerCell[fid];
			dint cell_n	= mesh->faceToNeighborCell[fid];

			dfloat fxUP = mymax(sgn(Flux[fid]),0.0);
			dfloat fxCD = mesh->fx[fid];

			dfloat lambda  = 0.0;

			for (int dim = 0; dim < ndimension; dim++) {

				dint fid_dim   = dim + fid*ndimension;
				dint cell_dim  = dim + cell*ndimension;
				dint cell_ndim = dim + cell_n*ndimension;

				dfloat phiP =  phiCell[cell_dim];
				dfloat phiN =  phiCell[cell_ndim];

				dfloat UP =  fxUP*phiP + (1-fxUP)*phiN;
				dfloat CD =  fxCD*phiP + (1-fxCD)*phiN;

				phiFace[fid_dim] = (1.0-lambda)*UP + lambda*CD;

			}
		}
	}


	if(interpolationScheme == interfaceCompression){
		// Interface Compression Scheme - Suitable only for scalar quantity volume fraction
		for (dint face = 0; face < mesh->nInnerAndPeriodicFaces; face++){

			dint fid 	= mesh->innerAndPeriodicFaces[face];
			dint cell 	= mesh->faceToOwnerCell[fid];
			dint cell_n	= mesh->faceToNeighborCell[fid];

			dfloat fxUP = mymax(sgn(Flux[fid]),0.0);
			dfloat fxCD = mesh->fx[fid];

			for (int dim = 0; dim < ndimension; dim++) {

				dint fid_dim   = dim + fid*ndimension;
				dint cell_dim  = dim + cell*ndimension;
				dint cell_ndim = dim + cell_n*ndimension;

				dfloat phiP =  phiCell[cell_dim];
				dfloat phiN =  phiCell[cell_ndim];

				dfloat lambda  =   mymin(mymax(1 - mymax(sqr(1 - 4*phiP*(1 - phiP)), sqr(1 - 4*phiN*(1 - phiN))),0), 1); // limiter

				dfloat UP =  fxUP*phiP + (1-fxUP)*phiN;
				dfloat CD =  fxCD*phiP + (1-fxCD)*phiN;

				phiFace[fid_dim] = (1.0-lambda)*UP + lambda*CD;
			}

		}
	}


	if(interpolationScheme == Gamma){


		// This is Hard Coded for 2D
		if ((blendingCoeff < 0.0) || (blendingCoeff > 0.5)){
			printf ( "\n" );
			printf ( "Fatal error!\n" );
			printf ( " blendingCoeff should be >= 0 and <= 0.5" );
			exit ( 1 );
		}


		// GAMMA SCHEME

		dfloat *gradphi;		// gradient of field at cell centers
		createField(&gradphi,volumeVector,mesh);

		interpolateCellToFaceLinear(phiCell,phiFace,volumeScalar, mesh);

		fvdGradientOfScalar(1.0,phiFace,mesh,gradphi);




		for (dint face = 0; face < mesh->nInnerAndPeriodicFaces; face++){

			dint fid 	= mesh->innerAndPeriodicFaces[face];
			dint cell 	= mesh->faceToOwnerCell[fid];
			dint cell_n	= mesh->faceToNeighborCell[fid];

			dfloat gradphifdotd = 0.0;
			dfloat gradphiCdotd = 0.0;

			dfloat *gradCP;
			dfloat *gradCN;
			dfloat *d;

			gradCP = (dfloat*) calloc(mesh->ndimension, sizeof(dfloat));
			gradCN = (dfloat*) calloc(mesh->ndimension, sizeof(dfloat));
			d	   = (dfloat*) calloc(mesh->ndimension, sizeof(dfloat));


			for (dint dim = 0; dim < mesh->ndimension; dim++) {

				dint cell_dim  = dim + cell*ndimension;
				dint cell_ndim = dim + cell_n*ndimension;
				dint fid_dim   = dim + fid*ndimension;

				gradCP[dim] = gradphi[cell_dim];
				gradCN[dim] = gradphi[cell_ndim];
				d[dim] 	    = mesh->d[fid_dim];

			}

			dfloat phiP =  phiCell[cell];
			dfloat phiN =  phiCell[cell_n];

			gradphifdotd = phiN - phiP;


			if(Flux[fid] > 0){

				gradphiCdotd = dotProduct(d,gradCP);
			}
			else{
				gradphiCdotd = dotProduct(d,gradCN);
			}

			stabilise(gradphiCdotd, mesh);

			dfloat phiC = 1.0 - 0.5*gradphifdotd/gradphiCdotd;




			//===================================================!

			dfloat fxUP = mymax(sgn(Flux[fid]),0.0);
			dfloat fxCD = mesh->fx[fid];

			stabilise(blendingCoeff, mesh);

			dfloat lambda = mymin(mymax(phiC/blendingCoeff, 0), 1);

			dfloat UP =  fxUP*phiP + (1-fxUP)*phiN;
			dfloat CD =  fxCD*phiP + (1-fxCD)*phiN;

			phiFace[fid] = (1.0-lambda)*UP + lambda*CD;


			free(gradCP);
			free(gradCN);
			free(d);


		}


		free(gradphi);
	}








	// Update Boundary Values from Boundary Cells

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
