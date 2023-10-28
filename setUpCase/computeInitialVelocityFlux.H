#include "mesh.h"

//#include "../setUpCase/computeInitialVelocityFlux.H"

void computeInitialVelocityFlux(mesh_t *mesh);

void computeInitialVelocityFlux(mesh_t *mesh){

	dfloat *UFace;
	createField(&UFace,surfaceVector,mesh);

	interpolateCellToFaceLinear(mesh->U,UFace,volumeVector,mesh);
	constructFluxofSurfaceVector(UFace,mesh->phiU,mesh);


	free(UFace);
}

