#include "mesh.h"

//#include "../Library/fvs/updateMixtureViscosity.H"

void updateMixtureViscosity(mesh_t *mesh);

void updateMixtureViscosity(mesh_t *mesh){

	if(mesh->alphaL > 0.0){ // Case for Liquid volume fraction is 1

		for (int cell = 0; cell < mesh->ncells; cell++) {
			mesh->mu[cell]  =  mesh->alpha[cell]*mesh->muL + (1-mesh->alpha[cell])*mesh->muG;
		}

	} else { // Case for Gas volume fraction is 1

		for (int cell = 0; cell < mesh->ncells; cell++) {
			mesh->mu[cell]  =  mesh->alpha[cell]*mesh->muG + (1-mesh->alpha[cell])*mesh->muL;
		}
	}

	computeDynamicViscosityFaceValues(mesh); //OK

}
