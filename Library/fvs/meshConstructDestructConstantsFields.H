#include "mesh.h"
//#include "../Library/fvs/meshConstructDestructConstantsFields.H"

void meshConstructConstantsFields(mesh_t *mesh);
void meshDestructConstantsFields(mesh_t *mesh);

void meshConstructConstantsFields(mesh_t *mesh){

	createField(&mesh->g,volumeVector,mesh);
	createField(&mesh->gf,surfaceVector,mesh);

}

void meshDestructConstantsFields(mesh_t *mesh){

	free(mesh->g);
	free(mesh->gf);


}
