#include "mesh.h"

//#include "../Library/fvs/updateHydroStaticPressure.H"

void updateHydroStaticPressure(mesh_t *mesh);

void updateHydroStaticPressure(mesh_t *mesh){


	if (mesh->includeGravity == false) {

		return;
	}



	// Estimate Hydro-Static pressure

	//===========================================================================!

	// Computing hydrostatic pressure

	for(int cell = 0;cell < mesh->ncells; ++cell){

		mesh->ph[cell] =  mesh->rho[cell]*mesh->gDotdh[cell];

	}



	// for North boundary Cells -- No integration

	// Internal Cells

	for( dint k = (mesh->ncellsz-1); k >=0; k--){
		for( dint i = 0; i < mesh->ncellsx; i++){
			for( dint j = 0; j < mesh->ncellsy; j++){

				dint cell = computeLinearIndex(i,j,k, mesh->ncellsx, mesh->ncellsy,mesh->ncellsz);

				dint nface  = mesh->topFace ; // north face

				dint fid  = nface + cell*mesh->ncellFaces;
				dint cell_n	=  mesh->faceToNeighborCell[fid];

				mesh->ph[cell] = mesh->ph[cell] + mesh->ph[cell_n];

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

			mesh->ph[bcell] = mesh->ph[bcell] + mesh->ph[ocell];



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

			mesh->ph[bcell] = mesh->ph[ocell];

		}
	}

	// for East boundary Cells


	for( dint j = 0; j < mesh->ncellsy; j++){
		for( dint k = 0; k < mesh->ncellsz; k++){

			dint i = (mesh->ncellsx-1);

			dint cell = computeLinearIndex(i,j,k, mesh->ncellsx,mesh->ncellsy,mesh->ncellsz);

			dint face = mesh->frontFace;
			dint fid  = face + cell*mesh->ncellFaces;

			int ocell =  mesh->faceToOwnerCell[fid];
			int bcell =  mesh->faceToBoundaryCell[fid];

			mesh->ph[bcell] = mesh->ph[ocell];

		}
	}


	// ZeroGradient for Axisymmetric Boundaries

	correctExtrapolationToAxisymmetricBoundaryCellScalar(mesh->ph,mesh);


//	strcpy(mesh->variableName, "ph");
//	printVolumeScalar(mesh->ph, mesh);


}
