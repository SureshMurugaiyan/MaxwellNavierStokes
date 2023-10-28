#include "mesh.h"
#include "../setUpCase/BoundaryNames.H"

void meshCalculateInnerFaceIDIncludingPeriodicFaces(mesh_t *mesh);

void meshCalculateInnerFaceIDIncludingPeriodicFaces(mesh_t *mesh){

	free(mesh->boundaryFaces);

	mesh->boundaryFaces = (dint*) calloc(mesh->nBoundaryFaces, sizeof(dint));

	for (dint i = 0; i < mesh->nBoundaryFaces; i++) {
		mesh->boundaryFaces[i] = 0;
	}


	dint i = 0;
	dint j = 0;

	for (int cell = 0; cell < mesh->nicells; cell++) {
		for (int face = 0; face < mesh->ncellFaces; face++) {

			int fid = face+cell*mesh->ncellFaces;
			int cell_n = mesh->faceToNeighborCell[fid];

			if (cell_n < mesh->nicells) {

				mesh->innerAndPeriodicFaces[i] = fid;

				i++;

			} else{

				mesh->boundaryFaces[j] = fid;

				j++;

			}


		}
	}



//	printf ( "!--------------------------------------------- !\n" );
//	printf ( "! Boundary Face ID                             !\n" );
//	printf ( "!--------------------------------------------- !\n" );
//
//	for (dint face = 0; face < mesh->nBoundaryFaces; face++) {
//
//		printf ( " %5d \n", mesh->boundaryFaces[face] );
//
//	}
//
//
//	printf ( "!--------------------------------------------- !\n" );
//	printf ( "! Inner and Periodic Face ID                   !\n" );
//	printf ( "!--------------------------------------------- !\n" );
//
//	for (int face = 0; face < mesh->nInnerAndPeriodicFaces; face++) {
//
//		printf ( " %5d %5d \n", face , mesh->innerAndPeriodicFaces[face] );
//
//	}
//	printf ( "!--------------------------------------------- !\n" );

}
