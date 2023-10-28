#include "mesh.h"

void meshCalculateNoOfBoundaryandInnerFaces(mesh_t *mesh);

void meshCalculateNoOfBoundaryandInnerFaces(mesh_t *mesh){

	mesh->nBoundaryFaces = 0;
	mesh->nInnerFaces	 = 0;

	for (dint cell = 0; cell < mesh->nicells; cell++) {
		for (dint face = 0; face < mesh->ncellFaces; face++) {

			dint fid = face + cell*mesh->ncellFaces;

			if ( mesh->cellToFaces[fid]==(-1)){

				mesh->nBoundaryFaces = mesh->nBoundaryFaces+1;
			}

		}
	}


	mesh->nInnerFaces = mesh->nicells*mesh->ncellFaces-mesh->nBoundaryFaces;



//	printf ( "!--------------------------------------------- !\n" );
//	printf ( "! Number of Boundary and Inner Faces           !\n" );
//	printf ( "!--------------------------------------------- !\n!" );
//
//	printf(" No of Boundary Faces,     N = %d\t \n!", mesh->nBoundaryFaces);
//
//	printf(" No of Inner Faces,        N = %d\t \n!",mesh->nInnerFaces );
//	printf ( "--------------------------------------------- !\n" );

}
