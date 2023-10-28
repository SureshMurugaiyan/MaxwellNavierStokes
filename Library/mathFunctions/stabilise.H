#include "mesh.h"

dfloat stabilise(dfloat a, mesh_t *mesh);

dfloat stabilise(dfloat a, mesh_t *mesh){

	if(a < 0){
		a =  a - mesh->SMALL;
	}
	else{
		a =  a + mesh->SMALL;
	}

	return a;
}
