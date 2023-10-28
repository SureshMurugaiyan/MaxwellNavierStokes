#include "mesh.h"

//#include "../setUpCase/setConstants.H"

void setConstants(mesh_t *mesh);

void setConstants(mesh_t *mesh){

	mesh->rhoG 	 = 1.0;
	mesh->rhoL 	 = 1.0;

	mesh->muG 	 = 1.0;
	mesh->muL 	 = 1.0;

	mesh->Re  	 = 100.0;
	mesh->Fr  	 = 1.0;
	mesh->Ha  	 = 1.0;

	mesh->alphaG =  0.0;
	mesh->alphaL =  1.0;

	mesh->sigmaG = 1.0;
	mesh->sigmaL = 1.0;

	// Cell center values of g

	for(int cell = 0; cell < mesh->ncells; ++cell){

		dint idx	= mesh->dim_x + cell*mesh->ndimension;
		dint idy	= mesh->dim_y + cell*mesh->ndimension;
		dint idz	= mesh->dim_z + cell*mesh->ndimension;

		mesh->g[idx]	=    0.0;
		mesh->g[idy]	=    0.0;
		mesh->g[idz]	=   -1.0;

	}

	// Cell face values of g

	for (int cell = 0; cell < mesh->nicells; cell++) {

		for (int face = 0; face < mesh->ncellFaces; face++) {

			dint fid 	= face + cell*mesh->ncellFaces;

			dint fidx	=  mesh->dim_x + fid*mesh->ndimension;
			dint fidy	=  mesh->dim_y + fid*mesh->ndimension;
			dint fidz	=  mesh->dim_z + fid*mesh->ndimension;

			mesh->gf[fidx]	=   0.0;
			mesh->gf[fidy]	=   0.0;
			mesh->gf[fidz]	=  -1.0;

		}

	}


//	strcpy(mesh->variableName, "g");
//	printVolumeVector(mesh->g, mesh);
//
//	strcpy(mesh->variableName, "gf");
//	printSurfaceVector(mesh->gf, mesh);

}

