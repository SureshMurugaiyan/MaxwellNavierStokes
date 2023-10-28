#include "mesh.h"

//#include "../Library/fvs/storeReferencePressure.H"

void storeReferencePressure(mesh_t *mesh);

void storeReferencePressure(mesh_t *mesh){

	for (dint cell = 0; cell < mesh->ncells; cell++) {

		mesh->pstar[cell]   = mesh->p[cell];

	}

}
