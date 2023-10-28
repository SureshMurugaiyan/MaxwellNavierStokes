#include "mesh.h"

//#include "../Library/fvc/fvcDivergenceOfScalar.H"

dfloat fvcDivergenceOfScalar
(
	dfloat GD,
	dfloat* Phif,
	dfloat* Flux,
	mesh_t *mesh
);

dfloat fvcDivergenceOfScalar
(
	dfloat GD,
	dfloat* Phif,
	dfloat* Flux,
	mesh_t *mesh
)
{

	for (int cell = 0; cell < mesh->nicells; cell++) {
		for (int face = 0; face < mesh->ncellFaces; face++) {

			dint fid 	   = face + cell*mesh->ncellFaces;

			dfloat R_new   = - GD*mesh->rV[cell]*Phif[fid]*Flux[fid];

			mesh->bP[cell] =	mesh->bP[cell] + R_new;

		}
	}

	return 0.0;

}
