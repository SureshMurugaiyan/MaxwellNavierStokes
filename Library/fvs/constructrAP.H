#include "mesh.h"

//#include "../Library/fvs/constructrAP.H"

void constructrAP(mesh_t *mesh);

void constructrAP(mesh_t *mesh){

	//  computeRAPCellandFaceValue

	for (int cell = 0; cell < mesh->nicells; cell++) {
		for (int dim = 0; dim < mesh->ndimension; dim++) {

			int cell_dim = dim + cell*mesh->ndimension;

			mesh->rAP[cell_dim] = 1.0/mesh->AAP[cell_dim];

		}
	}

	// extrapolate to boundary cells

	extrapolateToBoundaryCellVector(mesh->rAP,mesh);

	// ZeroGradient for Axisymmetric Boundaries

	correctExtrapolationToAxisymmetricBoundaryCellVector(mesh->rAP,mesh);


	// Computing cell face  values
	interpolateCellToFaceLinear(mesh->rAP,mesh->rAPf,volumeVector,mesh);

}


