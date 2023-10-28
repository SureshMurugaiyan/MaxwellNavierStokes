#include "mesh.h"

//#include "../Library/fvs/constructDestructSolutionFields.H"

void constructSolutionFields(mesh_t *mesh);
void destructSolutionFields(mesh_t *mesh);

void constructSolutionFields(mesh_t *mesh){

	createField(&mesh->U,volumeVector3,mesh);
	createField(&mesh->p,volumeScalar3,mesh);
	createField(&mesh->alpha,volumeScalar3,mesh);

	createField(&mesh->rho,volumeScalar,mesh);
	createField(&mesh->mu,volumeScalar,mesh);

	createField(&mesh->rhof,surfaceScalar,mesh);
	createField(&mesh->muf,surfaceScalar,mesh);

	createField(&mesh->phiU,surfaceScalar,mesh);
	createField(&mesh->rhophiU,surfaceScalar,mesh);


	createField(&mesh->Cof,surfaceScalar,mesh);
	createField(&mesh->Co,volumeScalar,mesh);

	createField(&mesh->Ustar,volumeVector,mesh);
	createField(&mesh->pstar,volumeScalar,mesh);

	createField(&mesh->rAP,volumeVector,mesh);
	createField(&mesh->rAPf,surfaceVector,mesh);

	createField(&mesh->HbyA,volumeVector,mesh);
	createField(&mesh->phiHbyA,surfaceScalar,mesh);

	createField(&mesh->ph,volumeScalar,mesh);

	//Pressure gradient variables
	createField(&mesh->gradPbyA,volumeVector,mesh);
	createField(&mesh->phigradPbyA,surfaceScalar,mesh);

	//Body force variables
	createField(&mesh->BF1,volumeVector,mesh);
	createField(&mesh->BF2,volumeVector,mesh);
	createField(&mesh->BF3,volumeVector,mesh);

	createField(&mesh->phi,volumeScalar,mesh);

	createField(&mesh->B,volumeVector,mesh);
	createField(&mesh->Bf,surfaceVector,mesh);
	createField(&mesh->phiB,surfaceScalar,mesh);

	createField(&mesh->sigma,volumeScalar,mesh);
	createField(&mesh->sigmaf,surfaceScalar,mesh);

	createField(&mesh->sigmaUcrossBFlux,surfaceScalar,mesh);
	createField(&mesh->sigmagradPhiFlux,surfaceScalar,mesh);
	createField(&mesh->phiJ,surfaceScalar,mesh);

	createField(&mesh->Jphi,volumeVector,mesh);
	createField(&mesh->Ju,volumeVector,mesh);
	createField(&mesh->J,volumeVector,mesh);

	createField(&mesh->A,volumeVector,mesh);
	createField(&mesh->B_J,volumeVector,mesh);

}

void destructSolutionFields(mesh_t *mesh){

	free(mesh->U);
	free(mesh->p);
	free(mesh->alpha);

	free(mesh->rho);
	free(mesh->mu);

	free(mesh->rhof);
	free(mesh->muf);

	free(mesh->phiU);
	free(mesh->rhophiU);

	free(mesh->Cof);
	free(mesh->Co);

	free(mesh->Ustar);
	free(mesh->pstar);

	free(mesh->rAP);
	free(mesh->rAPf);

	free(mesh->HbyA);
	free(mesh->phiHbyA);

	free(mesh->ph);

	//Pressure gradient variables
	free(mesh->gradPbyA);
	free(mesh->phigradPbyA);

	//Body force variables
	free(mesh->BF1);
	free(mesh->BF2);
	free(mesh->BF3);

	free(mesh->phi);

	free(mesh->B);
	free(mesh->Bf);
	free(mesh->phiB);

	free(mesh->sigma);
	free(mesh->sigmaf);

	free(mesh->sigmaUcrossBFlux);
	free(mesh->sigmagradPhiFlux);
	free(mesh->phiJ);

	free(mesh->Jphi);
	free(mesh->Ju);
	free(mesh->J);

	free(mesh->A);
	free(mesh->B_J);

}
