#include "mesh.h"

void createField(dfloat** phi, dint variableType, mesh_t *mesh);

void createField(dfloat** phi, dint variableType, mesh_t *mesh){

	dint n=0;

	switch (variableType)
	{

	case vertexScalar:	n = mesh->nvertices;break;
	case vertexVector:	n = mesh->nvertices*mesh->ndimension;break;

	}



	*phi = (dfloat*) calloc(n, sizeof(dfloat));

	for (int i = 0; i < n; i++) {
		(*phi)[i] = 0.0;
	}


}
