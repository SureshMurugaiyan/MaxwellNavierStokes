#include "mesh.h"

void meshCellToVerticesConnectivity(mesh_t *mesh);

void meshCellToVerticesConnectivity(mesh_t *mesh){

	int *globalvertices;
	globalvertices	= (int*) calloc(mesh->ncellVertices, sizeof(int));

	for(dint i = 0; i < mesh->ncellsx; i++){
		for(dint j = 0;j < mesh->ncellsy; j++){
			for(dint k = 0;k < mesh->ncellsz; k++){

				dint cellid = computeLinearIndex(i,j,k,mesh->ncellsx, mesh->ncellsy, mesh->ncellsz);

				globalvertices[0] =  computeLinearIndex(i,j,k,mesh->nverticesx, mesh->nverticesy, mesh->nverticesz);
				globalvertices[1] =  computeLinearIndex(i,j,k+1,mesh->nverticesx, mesh->nverticesy, mesh->nverticesz);
				globalvertices[2] =  computeLinearIndex(i,j+1,k,mesh->nverticesx, mesh->nverticesy, mesh->nverticesz);
				globalvertices[3] =  computeLinearIndex(i,j+1,k+1,mesh->nverticesx, mesh->nverticesy, mesh->nverticesz);
				globalvertices[4] =  computeLinearIndex(i+1,j,k,mesh->nverticesx, mesh->nverticesy, mesh->nverticesz);
				globalvertices[5] =  computeLinearIndex(i+1,j,k+1,mesh->nverticesx, mesh->nverticesy, mesh->nverticesz);
				globalvertices[6] =  computeLinearIndex(i+1,j+1,k,mesh->nverticesx, mesh->nverticesy, mesh->nverticesz);
				globalvertices[7] =  computeLinearIndex(i+1,j+1,k+1,mesh->nverticesx, mesh->nverticesy, mesh->nverticesz);


				for(int vert = 0;vert < mesh->ncellVertices; vert++){

					int vid = vert + cellid*mesh->ncellVertices;

					mesh->cellToVertex[vid] = globalvertices[vert];

				}

			}
		}
	}

	free(globalvertices);



//	printf ( "!--------------------------------------------- !\n" );
//	printf ( "! Cell To Vertex Connectivity                  !\n" );
//	printf ( "!--------------------------------------------- !\n" );
//
//	printf ( "%11s%6s%6s%6s%6s\n","! Vert ID :","NWest","SWest","SEast","NEast");
//	printf ( "!--------------------------------------------- !\n!" );
//	for (int cell = 0; cell < mesh->ncells; cell++) {
//
//		printf ( " Cell%3d :", cell );
//
//		for (int vert = 0; vert < mesh->ncellVertices; vert++) {
//
//			int vid = vert 	+ cell*mesh->ncellVertices;
//
//			printf ( "%6d", mesh->cellToVertex[vid] );
//		}
//
//		printf ( "\n!" );
//	}
//	printf ( "--------------------------------------------- !\n" );



}


