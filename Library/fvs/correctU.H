#include "mesh.h"
#include "../../setUpCase/FieldNames.H"

//#include "../Library/fvs/correctU.H"

void correctU(mesh_t *mesh);

void correctU(mesh_t *mesh){


	for (int cell = 0; cell < mesh->nicells; cell++) {
		for (int dim = 0; dim < mesh->ndimension; dim++) {

			int cell_dim = dim + cell*mesh->ndimension;

			mesh->U[cell_dim] = mesh->HbyA[cell_dim]- mesh->gradPbyA[cell_dim] ;


		}
	}


	updateBoundaryValues(mesh->U, U, mesh);



}


