#include "mesh.h"

void meshConstructFaceToCellConnectivityFields(mesh_t *mesh);
void meshDestructFaceToCellConnectivityFields(mesh_t *mesh);

void meshConstructFaceToCellConnectivityFields(mesh_t *mesh){


	mesh->innerFaces 	= (dint*) calloc(mesh->nInnerFaces, sizeof(dint));
	mesh->boundaryFaces = (dint*) calloc(mesh->nBoundaryFaces, sizeof(dint));

	for (dint i = 0; i < mesh->nInnerFaces; i++) {
		mesh->innerFaces[i] = 0;
	}

	for (dint i = 0; i < mesh->nBoundaryFaces; i++) {
		mesh->boundaryFaces[i] = 0;
	}


	createIntField(&mesh->faceToOwnerCell,surfaceScalar,mesh);
	createIntField(&mesh->faceToNeighborCell,surfaceScalar,mesh);



}

void meshDestructFaceToCellConnectivityFields(mesh_t *mesh){


	free(mesh->innerFaces);
	free(mesh->boundaryFaces);

	free(mesh->faceToOwnerCell);
	free(mesh->faceToNeighborCell);


}
