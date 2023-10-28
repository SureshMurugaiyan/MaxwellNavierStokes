#include "mesh.h"

void meshConstructVerticesFields(mesh_t *mesh);
void meshDestructVerticesFields(mesh_t *mesh);


void meshConstructVerticesFields(mesh_t *mesh){

	createField(&mesh->vertices,vertexVector,mesh);

	createCharField(&mesh->variableName,mesh);
	createCharField(&mesh->directoryName,mesh);
	createCharField(&mesh->subdirectoryName,mesh);
	createCharField(&mesh->fileName,mesh);

}

void meshDestructVerticesFields(mesh_t *mesh){

	free(mesh->vertices);

	free(mesh->variableName);
	free(mesh->directoryName);
	free(mesh->subdirectoryName);
	free(mesh->fileName);

}
