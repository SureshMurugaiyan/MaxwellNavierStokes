#include "mesh.h"

void meshConstructInnerAndPeriodicFacesFields(mesh_t *mesh);

void meshDestructInnerAndPeriodicFacesFields(mesh_t *mesh);

void meshConstructInnerAndPeriodicFacesFields(mesh_t *mesh){

	mesh->innerAndPeriodicFaces = (int*) calloc(mesh->nInnerAndPeriodicFaces, sizeof(int));

	for (int face = 0; face < mesh->nInnerAndPeriodicFaces; face++) {
		mesh->innerAndPeriodicFaces[face] = 0;
	}

}

void meshDestructInnerAndPeriodicFacesFields(mesh_t *mesh){

	free(mesh->innerAndPeriodicFaces);

}
