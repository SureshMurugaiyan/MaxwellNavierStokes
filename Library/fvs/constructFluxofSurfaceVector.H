#include "mesh.h"

//#include "../Library/fvs/constructFluxofSurfaceVector.H"

void constructFluxofSurfaceVector(dfloat* PhiFace, dfloat* Flux, mesh_t *mesh);

void constructFluxofSurfaceVector(dfloat* PhiFace, dfloat* Flux, mesh_t *mesh){

	for (int cell = 0; cell < mesh->nicells; cell++) {
		for (int face = 0; face < mesh->ncellFaces; face++) {

			dint fid = face+cell*mesh->ncellFaces;

			dfloat sum = 0.0;

			for (int dim = 0; dim < mesh->ndimension; dim++) {

				int fid_dim = dim+ fid*mesh->ndimension;

				sum	=	sum + PhiFace[fid_dim]*mesh->Sf[fid_dim];

			}

			Flux[fid] = sum;
		}
	}

}
