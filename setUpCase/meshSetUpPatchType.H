
#include "mesh.h"
#include "BoundaryNames.H"

void meshSetUpPatchType(mesh_t *mesh);

void meshSetUpPatchType(mesh_t *mesh){

	// choose one of the 3 options
	// 1. emptyPatch
	// 2. periodicPatch
	// 3. periodicPatch

	for (int patch = 0; patch < mesh->npatchs; patch++) {

		switch (mesh->patchID[patch])
		{
		case topBoundary:
			mesh->patchType[patch]			=	boundaryPatch;
			break;
		case bottomBoundary:
			mesh->patchType[patch]			=	boundaryPatch;
			break;
		case leftBoundary:
			mesh->patchType[patch]			=	boundaryPatch;
			break;
		case rightBoundary:
			mesh->patchType[patch]			=	boundaryPatch;
			break;
		case frontBoundary:
			mesh->patchType[patch]			=	boundaryPatch;
			break;
		case backBoundary:
			mesh->patchType[patch]			=	boundaryPatch;
			break;
		}

	}

//	for (int patch = 0; patch < mesh->npatchs; patch++) {
//
//		switch (mesh->patchID[patch])
//		{
//		case topBoundary:
//			mesh->patchType[patch]			=	periodicPatch;
//			break;
//		case bottomBoundary:
//			mesh->patchType[patch]			=	periodicPatch;
//			break;
//		case leftBoundary:
//			mesh->patchType[patch]			=	periodicPatch;
//			break;
//		case rightBoundary:
//			mesh->patchType[patch]			=	periodicPatch;
//			break;
//		case frontBoundary:
//			mesh->patchType[patch]			=	periodicPatch;
//			break;
//		case backBoundary:
//			mesh->patchType[patch]			=	periodicPatch;
//			break;
//		}
//
//	}

//
//	printf ( "!--------------------------------------------- !\n" );
//	printf ( "! patchType                                    !\n" );
//	printf ( "!--------------------------------------------- !\n!" );
//
//	printf ( " %16s  %8s \n!","Patch ID","patchType" );
//	printf ( "--------------------------------------------- !\n!" );
//	for (dint patch = 0; patch < mesh->npatchs; patch++) {
//		printf ( " %18s  %18s \n!", getPatchID(mesh->patchID[patch]) ,
//				getPatchType(mesh->patchType[patch]) );
//	}
//
//	printf ( "--------------------------------------------- !\n" );

}
