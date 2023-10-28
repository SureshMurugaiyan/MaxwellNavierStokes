#include "mesh.h"

void assignMeshDimensions(mesh_t *mesh);

void assignMeshDimensions(mesh_t *mesh){

	mesh->ncells  = mesh->ncellsx*mesh->ncellsy*mesh->ncellsz;

	mesh->nverticesx = mesh->ncellsx + 1;
	mesh->nverticesy = mesh->ncellsy + 1;
	mesh->nverticesz = mesh->ncellsz + 1;

	mesh->nvertices  = mesh->nverticesx*mesh->nverticesy*mesh->nverticesz;

	mesh->ndimension 	= 3;

	mesh->ncellVertices = 8;


//	printf ( "!--------------------------------- !\n" );
//	printf ( "! assignMeshDimensions             !\n" );
//	printf ( "!--------------------------------- !\n!" );
//	printf ( "%16s %5d \n!","ncellsx = ",mesh->ncellsx);
//	printf ( "%16s %5d \n!","ncellsy = ",mesh->ncellsy);
//	printf ( "%16s %5d \n!","ncellsz = ",mesh->ncellsz);
//	printf ( "%16s %5d \n!","ncells = ",mesh->ncells);
//	printf ( "\n!");
//
//	printf ( "%16s %5d \n!","nverticesx = ",mesh->nverticesx);
//	printf ( "%16s %5d \n!","nverticesy = ",mesh->nverticesy);
//	printf ( "%16s %5d \n!","nverticesz = ",mesh->nverticesz);
//	printf ( "%16s %5d \n!","nvertices = ",mesh->nvertices);
//	printf ( "\n!");
//	printf ( "%16s %5d \n!","ndimension = ",mesh->ndimension);
//	printf ( "%16s %5d \n","ncellVertices = ",mesh->ncellVertices);
//	printf ( "!--------------------------------- !\n" );



}


