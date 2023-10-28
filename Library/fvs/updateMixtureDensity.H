#include "mesh.h"

//#include "../Library/fvs/updateMixtureDensity.H"

void updateMixtureDensity(mesh_t *mesh);

void updateMixtureDensity(mesh_t *mesh){

	if(mesh->alphaL > 0.0){ // Case for Liquid volume fraction is 1

		for (int cell = 0; cell < mesh->ncells; cell++) {
			mesh->rho[cell]  =  mesh->alpha[cell]*mesh->rhoL + (1-mesh->alpha[cell])*mesh->rhoG;
		}

	} else { // Case for Gas volume fraction is 1

		for (int cell = 0; cell < mesh->ncells; cell++) {
			mesh->rho[cell]  =  mesh->alpha[cell]*mesh->rhoG + (1-mesh->alpha[cell])*mesh->rhoL;
		}
	}

	computeDensityFaceValues(mesh); //OK

}
