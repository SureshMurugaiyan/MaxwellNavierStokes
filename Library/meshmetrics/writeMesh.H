#include "mesh.h"

void writeMesh(mesh_t *mesh);

void writeMesh(mesh_t *mesh){

	strcpy(mesh->variableName, "Vert");
	writeMeshData(mesh->vertices,vertexVector,mesh);

	strcpy(mesh->variableName, "Cf");
	writeMeshData(mesh->Cf,surfaceVector,mesh);

	strcpy(mesh->variableName, "Sf");
	writeMeshData(mesh->Sf,surfaceVector,mesh);

	strcpy(mesh->variableName, "magSf");
	writeMeshData(mesh->magSf,surfaceScalar,mesh);

	strcpy(mesh->variableName, "Sn");
	writeMeshData(mesh->Sn,surfaceVector,mesh);

	strcpy(mesh->variableName, "Vol");
	writeMeshData(mesh->V,volumeScalar,mesh);

	strcpy(mesh->variableName, "fx");
	writeMeshData(mesh->fx,surfaceScalar,mesh);

	strcpy(mesh->variableName, "d");
	writeMeshData(mesh->d,surfaceVector,mesh);

	strcpy(mesh->variableName, "magd");
	writeMeshData(mesh->magd,surfaceScalar,mesh);

	strcpy(mesh->variableName, "df");
	writeMeshData(mesh->df,surfaceVector,mesh);

	strcpy(mesh->variableName, "magdf");
	writeMeshData(mesh->magdf,surfaceScalar,mesh);

	strcpy(mesh->variableName, "delta");
	writeMeshData(mesh->delta,surfaceScalar,mesh);

}
