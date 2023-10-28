#include "mesh.h"

void meshConstructComputationalDomainFields(mesh_t *mesh);
void meshDestructComputationalDomainFields(mesh_t *mesh);

void meshConstructComputationalDomainFields(mesh_t *mesh){

	mesh->xi = (dfloat*) calloc( mesh->nverticesx, sizeof(dfloat));

	for (dint i = 0; i < mesh->nverticesx; i++) {
		(mesh->xi)[i] = 0.0;
	}


	mesh->eta = (dfloat*) calloc( mesh->nverticesy, sizeof(dfloat));

	for (dint i = 0; i < mesh->nverticesy; i++) {
		(mesh->eta)[i] = 0.0;
	}

	mesh->zeta = (dfloat*) calloc( mesh->nverticesz, sizeof(dfloat));

	for (dint i = 0; i < mesh->nverticesz; i++) {
		(mesh->zeta)[i] = 0.0;
	}

}


void meshDestructComputationalDomainFields(mesh_t *mesh){

	free(mesh->xi);
	free(mesh->eta);
	free(mesh->zeta);

}
