#include "mesh.h"

//#include "../Library/fvs/computeExternalMagneticFieldFlux.H"

void computeExternalMagneticFieldFlux(mesh_t *mesh);

void computeExternalMagneticFieldFlux(mesh_t *mesh){

	interpolateCellToFaceLinear(mesh->B,mesh->Bf,volumeVector,mesh);

	constructFluxofSurfaceVector(mesh->Bf,mesh->phiB,mesh);

//	strcpy(mesh->variableName, "Bf");
//	printSurfaceVector(mesh->Bf, mesh);

//	strcpy(mesh->variableName, "phiB");
//	printSurfaceScalar(mesh->phiB, mesh);

}

