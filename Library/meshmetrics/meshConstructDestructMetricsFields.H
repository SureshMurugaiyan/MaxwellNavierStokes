#include "mesh.h"

void meshConstructMetricsFields(mesh_t *mesh);
void meshDestructMetricsFields(mesh_t *mesh);

void meshConstructMetricsFields(mesh_t *mesh){

	createField(&mesh->C,volumeVector,mesh);
	createField(&mesh->Cf,surfaceVector,mesh);
	createField(&mesh->Sf,surfaceVector,mesh);
	createField(&mesh->magSf,surfaceScalar,mesh);
	createField(&mesh->Sn,surfaceVector,mesh);

	createField(&mesh->V,volumeScalar,mesh);
	createField(&mesh->rV,volumeScalar,mesh);

	createField(&mesh->d,surfaceVector,mesh);
	createField(&mesh->magd,surfaceScalar,mesh);
	createField(&mesh->rmagd,surfaceScalar,mesh);

	createField(&mesh->df,surfaceVector,mesh);
	createField(&mesh->magdf,surfaceScalar,mesh);
	createField(&mesh->rmagdf,surfaceScalar,mesh);

	createField(&mesh->delta,surfaceScalar,mesh);
	createField(&mesh->rdelta,surfaceScalar,mesh);

	createField(&mesh->fx,surfaceScalar,mesh);

	mesh->boundaryFacePatchID
	= (dint*) calloc(mesh->nBoundaryFaces, sizeof(dint));

	for (dint i = 0; i < mesh->nBoundaryFaces; i++) {
		mesh->boundaryFacePatchID[i]     = 0;
	}


}

void meshDestructMetricsFields(mesh_t *mesh){

	free(mesh->C);
	free(mesh->Cf);
	free(mesh->Sf);
	free(mesh->magSf);
	free(mesh->Sn);

	free(mesh->V);
	free(mesh->rV);

	free(mesh->d);
	free(mesh->magd);
	free(mesh->rmagd);

	free(mesh->df);
	free(mesh->magdf);
	free(mesh->rmagdf);

	free(mesh->delta);
	free(mesh->rdelta);

	free(mesh->fx);

	free(mesh->boundaryFacePatchID);


}
