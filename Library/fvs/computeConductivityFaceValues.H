#include "mesh.h"

//#include "../Library/fvs/computeConductivityFaceValues.H"

void computeConductivityFaceValues(mesh_t *mesh);

void computeConductivityFaceValues(mesh_t *mesh){

	computeCellFaceValues(mesh->sigma, mesh->sigmaf,volumeScalar,mesh->phiU,mesh->divScheme,mesh->divSchemeblendFactor,mesh);

//	strcpy(mesh->variableName, "rhof");
//	printSurfaceScalar(mesh->rhof, mesh);

}


