#include "mesh.h"

//#include "../Library/fvs/updateMixtureConductivity.H"

void updateMixtureConductivity(mesh_t *mesh);

void updateMixtureConductivity(mesh_t *mesh){

	if(mesh->alphaL > 0.0){ // Case for Liquid volume fraction is 1

		for (int cell = 0; cell < mesh->ncells; cell++) {
			mesh->sigma[cell]  =  mesh->alpha[cell]*mesh->sigmaL + (1-mesh->alpha[cell])*mesh->sigmaG;
		}

	} else { // Case for Gas volume fraction is 1

		for (int cell = 0; cell < mesh->ncells; cell++) {
			mesh->sigma[cell]  =  mesh->alpha[cell]*mesh->sigmaG + (1-mesh->alpha[cell])*mesh->sigmaL;
		}
	}

	computeConductivityFaceValues(mesh); //OK

}
