#include "mesh.h"

//#include "../Library/fvd/fvdGradientOfScalar.H"

void fvdGradientOfScalar
(
	dfloat GG,
	dfloat *phif,
	mesh_t *mesh,
	dfloat* GradientOfScalar
);

void fvdGradientOfScalar
(
	dfloat GG,
	dfloat *phif,
	mesh_t *mesh,
	dfloat* GradientOfScalar
)
{

	dfloat *sum;
	sum	= (dfloat*) calloc(mesh->ndimension, sizeof(dfloat));


	for (int cell = 0; cell < mesh->nicells; cell++) {

		for (int i = 0; i < mesh->ndimension; i++) {sum[i]=0.0;}

		for (int face = 0; face < mesh->ncellFaces; face++) {
			for (int dim = 0; dim < mesh->ndimension; dim++) {

				dint fid 	 = face + cell*mesh->ncellFaces;
				dint fid_dim = dim + fid*mesh->ndimension;

				sum[dim]	 = sum[dim] + phif[fid]*mesh->Sf[fid_dim];

			}
		}

		for (int dim = 0; dim < mesh->ndimension; dim++) {

			int cell_dim = dim+ cell*mesh->ndimension;

			GradientOfScalar[cell_dim]  = GG*sum[dim]*mesh->rV[cell];

		}


	}

	free(sum);

}
