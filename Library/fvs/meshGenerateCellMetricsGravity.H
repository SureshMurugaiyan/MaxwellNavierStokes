#include "mesh.h"

//#include "../Library/fvs/meshGenerateCellMetricsGravity.H"

void meshGenerateCellMetricsGravity(mesh_t *mesh);

void meshGenerateCellMetricsGravity(mesh_t *mesh){


	// Estimate reference Co-ordinate X

	for (int cell = 0; cell < mesh->ncells; cell++) {
		for (dint dimension = 0; dimension < mesh->ndimension; dimension++) {

			dint cell_dim = dimension + cell*mesh->ndimension;

			mesh->X[cell_dim] =  mesh->C[cell_dim]-mesh->Xref[dimension];

		}
	}


	// Estimate reference Co-ordinate Xf

	// Cell face values of g

	for (int cell = 0; cell < mesh->nicells; cell++) {

		for (int face = 0; face < mesh->ncellFaces; face++) {

			dint fid 	= face + cell*mesh->ncellFaces;

			for (dint dimension = 0; dimension < mesh->ndimension; dimension++) {

				dint fid_dim = dimension + fid*mesh->ndimension;

				mesh->Xf[fid_dim] =  mesh->Cf[fid_dim]-mesh->Xref[dimension];

			}

		}

	}


	// Estimate g dot X

	for(int cell = 0;cell < mesh->ncells; ++cell){

		dint idx	= mesh->dim_x + cell*mesh->ndimension;
		dint idy	= mesh->dim_y + cell*mesh->ndimension;
		dint idz	= mesh->dim_z + cell*mesh->ndimension;

		dfloat x = mesh->X[idx];
		dfloat y = mesh->X[idy];
		dfloat z = mesh->X[idz];

		dfloat gx	= mesh->g[idx];
		dfloat gy	= mesh->g[idy];
		dfloat gz	= mesh->g[idz];

		mesh->gDotX[cell] = gx*x + gy*y + gz*z;
	}


	for (int cell = 0; cell < mesh->nicells; cell++) {

		for (int face = 0; face < mesh->ncellFaces; face++) {

			dint fid 	= face + cell*mesh->ncellFaces;

			dint fidx	=  mesh->dim_x + fid*mesh->ndimension;
			dint fidy	=  mesh->dim_y + fid*mesh->ndimension;
			dint fidz	=  mesh->dim_z + fid*mesh->ndimension;

			dfloat xf = mesh->Xf[fidx];
			dfloat yf = mesh->Xf[fidy];
			dfloat zf = mesh->Xf[fidz];

			dfloat gxf	= mesh->gf[fidx];
			dfloat gyf	= mesh->gf[fidy];
			dfloat gzf	= mesh->gf[fidz];

			mesh->gDotXf[fid]	= gxf*xf + gyf*yf + gzf*zf;


		}

	}

	// Computing del-h

	// Internal Cells

	for(int cell = 0;cell < mesh->nicells; ++cell){

		dint nface  = mesh->topFace ; // north face

		dint fid  = nface + cell*mesh->ncellFaces;
		dint cell_n	=  mesh->faceToNeighborCell[fid];

		for(dint dimension=0;dimension<mesh->ndimension;dimension++){

			dint cell_dim   = dimension + cell*mesh->ndimension;
			dint cell_n_dim = dimension + cell_n*mesh->ndimension;

			mesh->dh[cell_dim]	=  -(mesh->C[cell_n_dim]- mesh->C[cell_dim]);
		}
	}

	// for North boundary Cells

	for( dint i = 0; i < mesh->ncellsx; i++){
		for( dint j = 0; j < mesh->ncellsy; j++){

			dint k = mesh->ncellsz-1;

			dint cell = computeLinearIndex(i,j,k, mesh->ncellsx,mesh->ncellsy,mesh->ncellsz);

			dint face = mesh->topFace;
			dint fid  = face + cell*mesh->ncellFaces;

			int ocell =  mesh->faceToOwnerCell[fid];
			int bcell =  mesh->faceToBoundaryCell[fid];

			for(dint dimension=0;dimension<mesh->ndimension;dimension++){

				dint ocell_dim   = dimension + ocell*mesh->ndimension;
				dint bcell_dim   = dimension + bcell*mesh->ndimension;

//				mesh->dh[bcell_dim]	=  -(mesh->Xref[dimension] - mesh->C[bcell_dim]);

				mesh->dh[bcell_dim]	=  0.0;

			}

		}
	}




	// for south boundary Cells
	for( dint i = 0; i < mesh->ncellsx; i++){
		for( dint j = 0; j < mesh->ncellsy; j++){

			dint k = 0;

			dint cell = computeLinearIndex(i,j,k, mesh->ncellsx,mesh->ncellsy,mesh->ncellsz);

			dint face = mesh->bottomFace;
			dint fid  = face + cell*mesh->ncellFaces;

			int ocell =  mesh->faceToOwnerCell[fid];
			int bcell =  mesh->faceToBoundaryCell[fid];


			for(dint dimension=0;dimension<mesh->ndimension;dimension++){

				dint ocell_dim   = dimension + ocell*mesh->ndimension;
				dint bcell_dim   = dimension + bcell*mesh->ndimension;

				mesh->dh[bcell_dim]	=  -(mesh->C[ocell_dim]- mesh->C[bcell_dim]);
			}

		}
	}



	// for West boundary Cells


	for( dint j = 0; j < mesh->ncellsy; j++){
		for( dint k = 0; k < mesh->ncellsz; k++){

			dint i = 0;

			dint cell = computeLinearIndex(i,j,k, mesh->ncellsx,mesh->ncellsy,mesh->ncellsz);

			dint face = mesh->backFace;
			dint fid  = face + cell*mesh->ncellFaces;

			int ocell =  mesh->faceToOwnerCell[fid];
			int bcell =  mesh->faceToBoundaryCell[fid];

			for(dint dimension=0;dimension<mesh->ndimension;dimension++){

				dint ocell_dim   = dimension + ocell*mesh->ndimension;
				dint bcell_dim   = dimension + bcell*mesh->ndimension;

				mesh->dh[bcell_dim]	=  mesh->dh[ocell_dim];
			}

		}
	}


	//===================++COMPLETE++=====================


	// for East boundary Cells


	for( dint j = 0; j < mesh->ncellsy; j++){
		for( dint k = 0; k < mesh->ncellsz; k++){

			dint i = (mesh->ncellsx-1);

			dint cell = computeLinearIndex(i,j,k, mesh->ncellsx,mesh->ncellsy,mesh->ncellsz);

			dint face = mesh->frontFace;
			dint fid  = face + cell*mesh->ncellFaces;

			int ocell =  mesh->faceToOwnerCell[fid];
			int bcell =  mesh->faceToBoundaryCell[fid];

			for(dint dimension=0;dimension<mesh->ndimension;dimension++){

				dint ocell_dim   = dimension + ocell*mesh->ndimension;
				dint bcell_dim   = dimension + bcell*mesh->ndimension;

				mesh->dh[bcell_dim]	=  mesh->dh[ocell_dim];
			}

		}
	}


	// ZeroGradient for Axisymmetric Boundaries

	correctExtrapolationToAxisymmetricBoundaryCellVector(mesh->dh,mesh);





	// Compute g dot dh


	for(int cell = 0;cell < mesh->ncells; ++cell){

		dint idx	= mesh->dim_x + cell*mesh->ndimension;
		dint idy	= mesh->dim_y + cell*mesh->ndimension;
		dint idz	= mesh->dim_z + cell*mesh->ndimension;

		dfloat dhx = mesh->dh[idx];
		dfloat dhy = mesh->dh[idy];
		dfloat dhz = mesh->dh[idz];

		dfloat gx	= mesh->g[idx];
		dfloat gy	= mesh->g[idy];
		dfloat gz	= mesh->g[idz];

		mesh->gDotdh[cell] = gx*dhx + gy*dhy + gz*dhz;


	}














//	strcpy(mesh->variableName, "X");
//	printVolumeVector(mesh->X, mesh);

//	strcpy(mesh->variableName, "Xf");
//	printSurfaceVector(mesh->Xf, mesh);

//	strcpy(mesh->variableName, "gDotX");
//	printVolumeScalar(mesh->gDotX, mesh);

//	strcpy(mesh->variableName, "gDotXf");
//	printSurfaceScalar(mesh->gDotXf, mesh);

//	strcpy(mesh->variableName, "dh");
//	printVolumeVector(mesh->dh, mesh);

//	strcpy(mesh->variableName, "gDotdh");
//	printVolumeScalar(mesh->gDotdh, mesh);

}


