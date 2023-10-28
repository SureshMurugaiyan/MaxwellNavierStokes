#include "mesh.h"

//#include "../Library/fvc/fvcVolumetricSourceScalar.H"

dfloat fvcVolumetricSourceScalar
(
	dfloat GV,
	dfloat *phi,
	mesh_t *mesh
);

dfloat fvcVolumetricSourceScalar
(
	dfloat GV,
	dfloat *phi,
	mesh_t *mesh
)
{

	for (int cell = 0; cell < mesh->nicells; cell++) {

		dfloat r_new      =  -GV*phi[cell];

		mesh->bP[cell]  = mesh->bP[cell]  + r_new;

	}

	return 0.0;

}

