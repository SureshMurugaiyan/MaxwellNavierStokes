#include "mesh.h"

//#include "../Library/fvs/computeModifiedFaceArea.H"

void computeModifiedFaceArea(dfloat* alpha_f, dfloat* Sfm,dfloat* magSfm,mesh_t *mesh);

void computeModifiedFaceArea(dfloat* alpha_f, dfloat* Sfm,dfloat* magSfm,mesh_t *mesh){


	for( dint cell = 0; cell < mesh->nicells;cell++){
		for (dint face = 0; face < mesh->ncellFaces;face++) {

			dint fid = face + cell*mesh->ncellFaces;

			for(dint dimension = 0; dimension < mesh->ndimension; dimension++){

				dint fid_dim = dimension + fid*mesh->ndimension;

				Sfm[fid_dim] = mesh->Sf[fid_dim]*alpha_f[fid_dim];

			}

		}
	}


	dfloat *area;
	area = (dfloat*) calloc(mesh->ndimension, sizeof(dfloat));


	for( dint cell = 0; cell < mesh->nicells;cell++){
		for (dint face = 0; face < mesh->ncellFaces; ++face) {
			dint fid = face + cell*mesh->ncellFaces;

			for (dint dimension = 0; dimension < mesh->ndimension; dimension++) {
				dint fid_dim = dimension + fid*mesh->ndimension;
				area[dimension] = Sfm[fid_dim];
			}

			magSfm[fid] = sqrt(pow(area[0],2) + pow(area[1],2) + pow(area[2],2));
		}
	}

	free(area);

}
