#include "mesh.h"

//#include "../Library/fvs/computeDynamicViscosityFaceValues.H"

void computeDynamicViscosityFaceValues(mesh_t *mesh);

void computeDynamicViscosityFaceValues(mesh_t *mesh){

	interpolateCellToFaceLinear(mesh->mu,mesh->muf,volumeScalar,mesh);

//	strcpy(mesh->variableName, "muf");
//	printSurfaceScalar(mesh->muf, mesh);

}


