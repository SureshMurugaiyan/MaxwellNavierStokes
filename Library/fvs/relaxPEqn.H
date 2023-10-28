#include "mesh.h"

//#include "../Library/fvs/relaxPEqn.H"

void relaxPEqn(mesh_t *mesh);

void relaxPEqn(mesh_t *mesh){

	for (int cell = 0; cell < mesh->nicells; cell++) {


		mesh->bP[cell]  =  mesh->bP[cell] + (1-mesh->PEqnlambda)*(mesh->aP[cell]/mesh->PEqnlambda)*mesh->pstar[cell];


	}


	for (int cell = 0; cell < mesh->nicells; cell++) {


		mesh->aP[cell]  =  mesh->aP[cell] /mesh->PEqnlambda ;

	}

}
