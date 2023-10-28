#include "mesh.h"

void meshConstructDomainFields(mesh_t *mesh);
void meshDestructDomainFields(mesh_t *mesh);

void meshConstructDomainFields(mesh_t *mesh){

	createIntField(&mesh->domain,volumeScalar,mesh);
}

void meshDestructDomainFields(mesh_t *mesh){

	free(mesh->domain);

}
