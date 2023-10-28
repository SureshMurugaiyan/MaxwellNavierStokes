#include "mesh.h"

void createCharField(char** phi, mesh_t *mesh);

void createCharField(char** phi, mesh_t *mesh){

	dint n = 50;
	*phi = (char*) calloc(n, sizeof(char));

	for (dint i = 0; i < n; i++) {
		(*phi)[i] = '\0';
	}

}

















