#include "mesh.h"

#include "BoundaryNames.H"

void meshSetUpPeriodicNeighborPatches(mesh_t *mesh);

void meshSetUpPeriodicNeighborPatches(mesh_t *mesh){

	/*******************************************************************************************/
	/* This will have effect only if the patch type is set to periodic in "meshSetUpPatchType" */
	/*******************************************************************************************/

	for (int patch = 0; patch < mesh->npatchs; patch++) {

		switch (mesh->patchID[patch])
		{
		case topBoundary:
			mesh->periodicNeighborPatchID[patch]			=	bottomBoundary;
			break;
		case bottomBoundary:
			mesh->periodicNeighborPatchID[patch]			=	topBoundary;
			break;
		case leftBoundary:
			mesh->periodicNeighborPatchID[patch]			=	rightBoundary;
			break;
		case rightBoundary:
			mesh->periodicNeighborPatchID[patch]			=	leftBoundary;
			break;
		case frontBoundary:
			mesh->periodicNeighborPatchID[patch]			=	backBoundary;
			break;
		case backBoundary:
			mesh->periodicNeighborPatchID[patch]			=	frontBoundary;
			break;
		}

	}

//	printf ( "!--------------------------------------------- !\n" );
//	printf ( "! neighborPatch of Each Patch                  !\n" );
//	printf ( "!--------------------------------------------- !\n!" );
//
//	printf ( " %16s  %8s \n!","Patch ID","neighborPatch ID" );
//	printf ( "--------------------------------------------- !\n!" );
//	for (dint patch = 0; patch < mesh->npatchs; patch++) {
//		printf ( " %18s  %18s \n!", getPatchID(mesh->patchID[patch]) ,
//				getPatchID(mesh->periodicNeighborPatchID[patch]) );
//	}
//
//	printf ( "--------------------------------------------- !\n" );


}
