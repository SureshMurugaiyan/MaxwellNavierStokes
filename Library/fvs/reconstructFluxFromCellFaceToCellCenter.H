#include "mesh.h"

//#include "../Library/fvs/reconstructFluxFromCellFaceToCellCenter.H"

void reconstructFluxFromCellFaceToCellCenter(dfloat* Flux, dfloat* cellCenterVector,mesh_t *mesh);

void reconstructFluxFromCellFaceToCellCenter(dfloat* Flux, dfloat* cellCenterVector,mesh_t *mesh){

	dfloat *snFlux;
	dfloat *summedsnFlux;
	dfloat *sum;

	createField(&snFlux,surfaceVector,mesh);
	createField(&summedsnFlux,volumeVector,mesh);
	createField(&sum,vector,mesh);

	//============================================================================!
	for( dint cell = 0; cell < mesh->nicells;cell++){
		for (dint face = 0; face < mesh->ncellFaces; ++face) {
			for (int dim = 0; dim < mesh->ndimension; dim++) {

				dint fid    = face + cell*mesh->ncellFaces;
				dint fid_dim = dim+ fid*mesh->ndimension;

				snFlux[fid_dim] = mesh->Sn[fid_dim]*Flux[fid];

			}
		}
	}
	//============================================================================!

	for (dint cell = 0; cell < mesh->nicells; ++cell) {

		// Sum the face values
		for (dint n = 0; n < mesh->ndimension; ++n) {sum[n]=0;}

		for (dint face = 0; face < mesh->ncellFaces; ++face) {

			for (int dim = 0; dim < mesh->ndimension; dim++) {

				dint fid	=  face  + cell*mesh->ncellFaces;
				dint fid_dim = dim + fid*mesh->ndimension;

				sum[dim] = sum[dim] + snFlux[fid_dim];
			}
		}

		// Assign to cell center

		for (int dim = 0; dim < mesh->ndimension; dim++) {

			dint cell_dim  = dim + cell*mesh->ndimension;

			summedsnFlux[cell_dim] = sum[dim];
		}

	}


	//============================================================================!
	// Matrix Vector Product

	dint m,n;
	dint ma,mb,mc,md,me,mf,mg,mh,mi;
	dint va,vb,vc;

	for (dint cell = 0; cell < mesh->nicells; ++cell) {

		// matrix indices
		m=0;n=0;
		ma = n + m*mesh->ndimension + cell*mesh->ndimension*mesh->ndimension;

		m=0;n=1;
		mb = n + m*mesh->ndimension + cell*mesh->ndimension*mesh->ndimension;

		m=0;n=2;
		mc = n + m*mesh->ndimension + cell*mesh->ndimension*mesh->ndimension;

		m=1;n=0;
		md = n + m*mesh->ndimension + cell*mesh->ndimension*mesh->ndimension;

		m=1;n=1;
		me = n + m*mesh->ndimension + cell*mesh->ndimension*mesh->ndimension;

		m=1;n=2;
		mf = n + m*mesh->ndimension + cell*mesh->ndimension*mesh->ndimension;

		m=2;n=0;
		mg = n + m*mesh->ndimension + cell*mesh->ndimension*mesh->ndimension;

		m=2;n=1;
		mh = n + m*mesh->ndimension + cell*mesh->ndimension*mesh->ndimension;

		m=2;n=2;
		mi = n + m*mesh->ndimension + cell*mesh->ndimension*mesh->ndimension;

		// vector indices
		va = 0 + cell*mesh->ndimension;
		vb = 1 + cell*mesh->ndimension;
		vc = 2 + cell*mesh->ndimension;

		dfloat maa = mesh->inverseOfSummedSnDyadSf[ma];
		dfloat mbb = mesh->inverseOfSummedSnDyadSf[mb];
		dfloat mcc = mesh->inverseOfSummedSnDyadSf[mc];
		dfloat mdd = mesh->inverseOfSummedSnDyadSf[md];
		dfloat mee = mesh->inverseOfSummedSnDyadSf[me];
		dfloat mff = mesh->inverseOfSummedSnDyadSf[mf];
		dfloat mgg = mesh->inverseOfSummedSnDyadSf[mg];
		dfloat mhh = mesh->inverseOfSummedSnDyadSf[mh];
		dfloat mii = mesh->inverseOfSummedSnDyadSf[mi];

		dfloat vaa = summedsnFlux[va];
		dfloat vbb = summedsnFlux[vb];
		dfloat vcc = summedsnFlux[vc];


		cellCenterVector[va] = maa*vaa + mbb*vbb + mcc*vcc;
		cellCenterVector[vb] = mdd*vaa + mee*vbb + mff*vcc;
		cellCenterVector[vc] = mgg*vaa + mhh*vbb + mii*vcc;

	}

	free(snFlux);
	free(summedsnFlux);
	free(sum);


}

