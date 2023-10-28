#include "mesh.h"

//#include "../Library/fvs/computeCourantNumberUsingConvectiveVelocity.H"

void computeCourantNumberUsingConvectiveVelocity(mesh_t *mesh);

void computeCourantNumberUsingConvectiveVelocity(mesh_t *mesh){


	for (dint face = 0; face < mesh->nInnerAndPeriodicFaces; face++){

		dint fid 	= mesh->innerAndPeriodicFaces[face];

		dfloat *d;
		dfloat *Sf;

		d  = (dfloat*) calloc(mesh->ndimension, sizeof(dfloat));
		Sf = (dfloat*) calloc(mesh->ndimension, sizeof(dfloat));

		for (dint dim = 0; dim < mesh->ndimension; dim++) {
			dint fid_dim   = dim + fid*mesh->ndimension;
			d[dim] 	    = mesh->d[fid_dim];
			Sf[dim] 	= mesh->Sf[fid_dim];
		}

		mesh->Cof[fid] = magnitude(mesh->phiU[fid]*mesh->deltatime/stabilise(dotProduct(d,Sf),mesh));

		free(d);
		free(Sf);
	}

	for (dint face = 0; face < mesh->nBoundaryFaces; face++) {

		dint fid = mesh->boundaryFaces[face];

		dfloat *d;
		dfloat *Sf;

		d  = (dfloat*) calloc(mesh->ndimension, sizeof(dfloat));
		Sf = (dfloat*) calloc(mesh->ndimension, sizeof(dfloat));

		for (dint dim = 0; dim < mesh->ndimension; dim++) {
			dint fid_dim   = dim + fid*mesh->ndimension;
			d[dim] 	    = mesh->df[fid_dim];
			Sf[dim] 	= mesh->Sf[fid_dim];
		}

		mesh->Cof[fid] = magnitude(mesh->phiU[fid]*mesh->deltatime/stabilise(dotProduct(d,Sf),mesh));


		free(d);
		free(Sf);
	}


	// Maximum Cell Courant Number is maximum among 6 cell face Courant Number

	for( dint cell = 0; cell < mesh->nicells;cell++){

		dint leftFacefid  = mesh->leftFace + cell*mesh->ncellFaces;

		dfloat maxCourantNumber = mesh->Cof[leftFacefid];

		for (dint face = 1; face < mesh->ncellFaces; ++face) {

			dint fid  = face + cell*mesh->ncellFaces;

			maxCourantNumber = mymax(mesh->Cof[fid],maxCourantNumber);

		}

		mesh->Co[cell] = maxCourantNumber;
	}


	// Maximum  Courant Number in the domain

	mesh->CoU = mesh->Co[0];

	for( dint cell = 0; cell < mesh->nicells;cell++){

		mesh->CoU = mymax(mesh->Co[cell],mesh->CoU);

	}



}
