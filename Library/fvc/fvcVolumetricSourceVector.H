#include "mesh.h"

//#include "../Library/fvc/fvcVolumetricSourceVector.H"

dfloat fvcVolumetricSourceVector
(
	dfloat GV,
	dfloat *phi,
	mesh_t *mesh
);

dfloat fvcVolumetricSourceVector
(
	dfloat GV,
	dfloat *phi,
	mesh_t *mesh
)
{

	for (int cell = 0; cell < mesh->nicells; cell++) {
		for (dint dimension = 0; dimension < mesh->ndimension; dimension++) {

			dint cell_dim 		= dimension + cell*mesh->ndimension;

			dfloat r_new        =  -GV*phi[cell_dim];

			mesh->BP[cell_dim]  = mesh->BP[cell_dim]  + r_new;
		}
	}

	return 0.0;

}

