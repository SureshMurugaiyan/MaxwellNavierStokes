#include "mesh.h"

//#include "../Library/fvs/meshConstructDestructGravityMetricsFields.H"


void meshConstructGravityMetricsFields(mesh_t *mesh);
void meshDestructGravityMetricsFields(mesh_t *mesh);

void meshConstructGravityMetricsFields(mesh_t *mesh){

	createField(&mesh->Xref,vector,mesh);

	createField(&mesh->X,volumeVector,mesh);
	createField(&mesh->Xf,surfaceVector,mesh);

	createField(&mesh->gDotX,volumeScalar,mesh);
	createField(&mesh->gDotXf,surfaceScalar,mesh);

	createField(&mesh->dh,volumeVector,mesh);
	createField(&mesh->gDotdh,volumeScalar,mesh);

}

void meshDestructGravityMetricsFields(mesh_t *mesh){

	free(mesh->Xref);

	free(mesh->X);
	free(mesh->Xf);

	free(mesh->gDotX);
	free(mesh->gDotXf);

	free(mesh->dh);
	free(mesh->gDotdh);

}
