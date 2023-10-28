#include "mesh.h"

void meshConstructPatchConnectivityFields(mesh_t *mesh);
void meshDestructPatchConnectivityFields(mesh_t *mesh);

void meshConstructPatchConnectivityFields(mesh_t *mesh){

	mesh->patchID = (dint*) calloc(mesh->npatchs, sizeof(dint));
	for (int patch = 0; patch < mesh->npatchs; patch++) {
		mesh->patchID[patch] = 0;
	}

	mesh->npatchFaces = (dint*) calloc(mesh->npatchs, sizeof(int));
	for (int patch = 0; patch < mesh->npatchs; patch++) {
		mesh->npatchFaces[patch] = 0;
	}

	mesh->patchStartFaceIndex = (dint*) calloc(mesh->npatchs, sizeof(int));
	for (int patch = 0; patch < mesh->npatchs; patch++) {
		mesh->patchStartFaceIndex[patch] = 0;
	}


	mesh->patchToBoundaryFace = (dint*) calloc(mesh->nBoundaryFaces, sizeof(dint));
	for (dint i = 0; i < mesh->nBoundaryFaces; i++) {
		mesh->patchToBoundaryFace[i]     = 0;
	}

	mesh->patchType 	= (dint*) calloc(mesh->npatchs, sizeof(dint));

	for (dint patch = 0; patch < mesh->npatchs; patch++) {

		mesh->patchType[patch] = 0;
	}


	mesh->periodicNeighborPatchID 	= (dint*) calloc(mesh->npatchs, sizeof(dint));

	for (dint patch = 0; patch < mesh->npatchs; patch++) {

		mesh->periodicNeighborPatchID[patch] = 0;
	}

	mesh->periodicNeighborPatch	= (int*) calloc(mesh->npatchs, sizeof(int));

	for (int patch = 0; patch < mesh->npatchs; patch++) {
		mesh->periodicNeighborPatch[patch]= 0;
	}


	createIntField(&mesh->faceToBoundaryCell,surfaceScalar,mesh);


}

void meshDestructPatchConnectivityFields(mesh_t *mesh){

	free(mesh->patchID);
	free(mesh->npatchFaces);
	free(mesh->patchStartFaceIndex);
	free(mesh->patchToBoundaryFace);
	free(mesh->patchType);
	free(mesh->periodicNeighborPatchID);
	free(mesh->periodicNeighborPatch);
	free(mesh->faceToBoundaryCell);
}
