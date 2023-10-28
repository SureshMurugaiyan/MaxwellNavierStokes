#include "mesh.h"

#include "../setUpCase/FieldNames.H"

void constructAlphaEquation(mesh_t *mesh){


	reinitializeField(mesh->aP,volumeScalar,mesh);
	reinitializeField(mesh->aN,surfaceScalar,mesh);
	reinitializeField(mesh->bP,volumeScalar,mesh);

	dint alphaEqn;

	alphaEqn = fvmDdtOfScalar(1.0,mesh->alpha, alpha, euler,mesh);
	alphaEqn = fvmDivergenceOfScalar(1.0,mesh->alpha,alpha,mesh->phiU,mesh->divScheme,mesh->divSchemeblendFactor,mesh);


}

