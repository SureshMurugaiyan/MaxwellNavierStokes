#include "mesh.h"

void meshGenerateBoundaryPatchConnectivity(mesh_t *mesh);

void meshDeclareBoundaryPatchID(mesh_t *mesh);
void meshCountNoOfBoundaryFacesPerPatch(mesh_t *mesh);
void meshpatchStartFaceIndex(mesh_t *mesh);
void meshpatchToBoundaryFaceConnectivity(mesh_t *mesh);


void meshDeclareBoundaryPatchIDUnitTest(mesh_t *mesh);
void meshCountNoOfBoundaryFacesPerPatchUnitTest(mesh_t *mesh);
void meshpatchStartFaceIndexUnitTest(mesh_t *mesh);
void meshpatchToBoundaryFaceConnectivityUnitTest(mesh_t *mesh);

void meshGenerateBoundaryPatchConnectivity(mesh_t *mesh){

	meshDeclareBoundaryPatchID(mesh);

	meshCountNoOfBoundaryFacesPerPatch(mesh);

	meshpatchStartFaceIndex(mesh);

	meshpatchToBoundaryFaceConnectivity(mesh);

}



void meshDeclareBoundaryPatchID(mesh_t *mesh){

	int count = 0;
	int j = 0;

	// Use nested loop
	// Pick elements one by one
	for (int i = 0; i < mesh->nBoundaryFaces; i++) {
		// check if the element has appeared to the left side of array
		for (j = 0; j < i; j++){
			// if the elements has already appeared, break the for loop
			// and pick the next element
			if (mesh->boundaryFacePatchID[i] == mesh->boundaryFacePatchID[j])
			{
				break;
			}
		}

		// if the elements has already not appeared
		// then you will reach the condition i=j
		// increase the count by one
		if (i == j) {

			mesh->patchID[count] = mesh->boundaryFacePatchID[i];

			count++;


		}
	}


//	printf ( "!--------------------------------------------- !\n" );
//	printf ( "! patchID                                      !\n" );
//	printf ( "!--------------------------------------------- !\n" );
//
//	for (int patch = 0; patch < mesh->npatchs; patch++) {
//		printf ( "%4d %25s\n", mesh->patchID[patch],getPatchID(mesh->patchID[patch]) );
//		printf ( "\n");
//	}

}



void meshCountNoOfBoundaryFacesPerPatch(mesh_t *mesh){

	for (int patch = 0; patch < mesh->npatchs; patch++) {
		for (int face = 0; face < mesh->nBoundaryFaces; face++) {
			if (mesh->boundaryFacePatchID[face] == mesh->patchID[patch]){
				mesh->npatchFaces[patch] = mesh->npatchFaces[patch]+1;
			}
		}
	}


//	printf ( "!--------------------------------------------- !\n" );
//	printf ( "! Number of Boundary Faces In Each Patch       !\n" );
//	printf ( "!--------------------------------------------- !\n!" );
//
//	printf ( " %20s  %8s \n!","Patch ID","nFaces" );
//	printf ( "--------------------------------------------- !\n!" );
//	for (int patch = 0; patch < mesh->npatchs; patch++) {
//		printf ( " %20s  %8d \n!", getPatchID(mesh->patchID[patch]), mesh->npatchFaces[patch]);
//	}
//
//	printf ( "--------------------------------------------- !\n" );


}




void meshpatchStartFaceIndex(mesh_t *mesh){

	mesh->patchStartFaceIndex[0] = 0;

	int sum = 0;

	for (int patch = 1; patch < mesh->npatchs; patch++) {
		sum =  sum + mesh->npatchFaces[patch-1];

		mesh->patchStartFaceIndex[patch] = sum;
	}


//	printf ( "!--------------------------------------------------- !\n" );
//	printf ( "! Start Face Index In Each Patch                     !\n" );
//	printf ( "!--------------------------------------------------- !\n!" );
//
//	printf ( " %20s  %12s\n!","Patch ID","StartFaceIndex");
//	printf ( "--------------------------------------------------- !\n!" );
//	for (int patch = 0; patch < mesh->npatchs; patch++) {
//		printf ( " %20s  %8d  \n!", getPatchID(mesh->patchID[patch]),
//				mesh->patchStartFaceIndex[patch]);
//	}
//
//	printf ( "--------------------------------------------------- !\n" );


}


void meshpatchToBoundaryFaceConnectivity(mesh_t *mesh){

	int i = 0;
	for (int patch = 0; patch < mesh->npatchs; patch++) {
		for (int face = 0; face < mesh->nBoundaryFaces; face++) {
			if (mesh->boundaryFacePatchID[face] == mesh->patchID[patch]){
				mesh->patchToBoundaryFace[i]     = mesh->boundaryFaces[face];
				i=i+1;
			}
		}
	}



//	printf ( "!--------------------------------------------- !\n" );
//	printf ( "! Boundary Face ID in each Patch               !\n" );
//	printf ( "!--------------------------------------------- !\n" );
//
//
//	for (int patch = 0; patch < mesh->npatchs; patch++) {
//
//		printf ( "! Patch ID =  %4s \n!", getPatchID(mesh->patchID[patch]) );
//
//		int face = mesh->patchStartFaceIndex[patch];
//		printf ( " --------------------------------- \n!" );
//		for (int patchFace = 0; patchFace < mesh->npatchFaces[patch]; patchFace++) {
//
//			printf ( " Face ID =  %4d \n!", mesh->patchToBoundaryFace[face] );
//
//			face++;
//
//		}
//
//		printf ( " --------------------------------- \n" );
//	}
//
//
//	printf ( "!--------------------------------------------------- !\n" );
//	printf ( "! Start Face ID In Each Patch                        !\n" );
//	printf ( "!--------------------------------------------------- !\n!" );
//
//	printf ( " %20s  %12s\n!","Patch ID","StartFaceID");
//	printf ( "--------------------------------------------------- !\n!" );
//	for (int patch = 0; patch < mesh->npatchs; patch++) {
//		printf ( " %20s  %8d  \n!", getPatchID(mesh->patchID[patch]),
//				mesh->patchToBoundaryFace[mesh->patchStartFaceIndex[patch]]);
//	}
//
//	printf ( "--------------------------------------------------- !\n" );

}

