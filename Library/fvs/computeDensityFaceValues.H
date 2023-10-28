#include "mesh.h"

//#include "../Library/fvs/computeDensityFaceValues.H"

void computeDensityFaceValues(mesh_t *mesh);

void computeDensityFaceValues(mesh_t *mesh){

	computeCellFaceValues(mesh->rho, mesh->rhof,volumeScalar,mesh->phiU,mesh->divScheme,mesh->divSchemeblendFactor,mesh);

//	strcpy(mesh->variableName, "rhof");
//	printSurfaceScalar(mesh->rhof, mesh);

}


