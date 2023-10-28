#include "mesh.h"

void copyDataToFrom
(
	dfloat* copyTo,
	dfloat* copyFrom,
	int fieldType,
	mesh_t *mesh
);

void copyDataToFrom
(
	dfloat* copyTo,
	dfloat* copyFrom,
	int fieldType,
	mesh_t *mesh
)
{

	int n=0;

	switch (fieldType)
	{
	case vertexScalar:			n = mesh->ncells*mesh->ncellVertices;break;
	case vertexVector:			n = mesh->ncells*mesh->ncellVertices*mesh->ndimension;break;
	case volumeScalar:			n = mesh->ncells;break;
	case volumeVector:			n = mesh->ncells*mesh->ndimension;break;
	case surfaceScalar:			n = mesh->ncells*mesh->ncellFaces;break;
	case surfaceVector:			n = mesh->ncells*mesh->ncellFaces*mesh->ndimension;break;
	case scalar:				n = 1;break;
	case vector:				n = mesh->ndimension;break;
	case cellFaceVertexScalar:	n = mesh->ncellFaces*mesh->nfaceVertices;break;
	case cellFaceVertexVector:	n = mesh->ncellFaces*mesh->nfaceVertices*mesh->ndimension;break;
	}


	for (int i = 0; i< n; i++){
		copyTo[i]= copyFrom[i];
	}





//	printf("\n!");
//
//	printf ( "!--------------------------------------------- !\n" );
//	printf ( "! Copy Data                                    !\n" );
//	printf ( "!--------------------------------------------- !\n!");
//
//
//	for (int i = 0; i < n; i++) {
//
//		printf ( "%.3f,%.3f", copyFrom[i], copyTo[i]);
//		printf("\n!");
//
//	}
//
//	printf ( "--------------------------------------------- !\n" );

}

