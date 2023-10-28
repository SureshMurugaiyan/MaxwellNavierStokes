#include "mesh.h"

//#include "../Library/fvs/setReferencePressure.H"

void setReferencePressure(mesh_t *mesh);

void setReferencePressure(mesh_t *mesh){


	dfloat min_value = mesh->p[0];

	for (dint cell = 0; cell < mesh->ncells; cell++) {

		if( mesh->p[cell] < min_value) {
			min_value =  mesh->p[cell];
		}

	}

	if(min_value < 0) {

		for (dint cell = 0; cell < mesh->ncells; cell++) {
			mesh->p[cell] -= min_value;
		}

	}


}


