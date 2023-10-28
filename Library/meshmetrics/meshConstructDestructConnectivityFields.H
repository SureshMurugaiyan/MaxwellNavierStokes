#include "mesh.h"

void meshConstructConnectivityFields(mesh_t *mesh);
void meshDestructConnectivityFields(mesh_t *mesh);

void meshConstructConnectivityFields(mesh_t *mesh){

	createIntField(&mesh->cellToVertex,vertexScalar,mesh);
	createIntField(&mesh->cellToCells,surfaceScalar,mesh);
	createIntField(&mesh->cellToFaces,surfaceScalar,mesh);
	createIntField(&mesh->faceToVertices,cellFaceVertexScalar,mesh);

}

void meshDestructConnectivityFields(mesh_t *mesh){

	free(mesh->cellToVertex);
	free(mesh->cellToCells);
	free(mesh->cellToFaces);
	free(mesh->faceToVertices);

}
