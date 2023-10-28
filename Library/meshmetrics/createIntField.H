#include "mesh.h"

void createIntField(dint** phi, dint variableType, mesh_t *mesh);

void createIntField(dint** phi, dint variableType, mesh_t *mesh){

	dint n = 0;

	switch (variableType)
	{
	case vertexScalar:			n = mesh->ncells*mesh->ncellVertices;break;
	case vertexVector:			n = mesh->ncells*mesh->ncellVertices*mesh->ndimension;break;
	case volumeScalar:			n = mesh->ncells;break;
	case volumeScalar2:			n = mesh->ncells*2;break;
	case volumeScalar3:			n = mesh->ncells*3;break;
	case volumeVector:			n = mesh->ncells*mesh->ndimension;break;
	case volumeVector2:			n = mesh->ncells*mesh->ndimension*2;break;
	case volumeVector3:			n = mesh->ncells*mesh->ndimension*3;break;
	case surfaceScalar:			n = mesh->ncells*mesh->ncellFaces;break;
	case surfaceVector:			n = mesh->ncells*mesh->ncellFaces*mesh->ndimension;break;
	case scalar:			    n = 1;break;
	case vector:			    n = mesh->ndimension;break;
	case cellFaceVertexScalar:	n = mesh->ncellFaces*mesh->nfaceVertices;break;
	case cellFaceVertexVector:	n = mesh->ncellFaces*mesh->nfaceVertices*mesh->ndimension;break;
	case volumeTensor:			n = mesh->ncells*mesh->ndimension*mesh->ndimension;break;
	case surfaceTensor:			n = mesh->ncells*mesh->ncellFaces*mesh->ndimension*mesh->ndimension;break;
	}

	*phi = (dint*) calloc(n, sizeof(dint));

	for (dint i = 0; i < n; i++) {
		(*phi)[i] = 0;
	}


}
















