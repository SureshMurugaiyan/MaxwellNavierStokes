#include "mesh.h"

//#include "../Library/fvs/writeResults.H"

void writeResults(mesh_t *mesh);

void writeResults(mesh_t *mesh){


	//=============================================================================!
	// Primary Fields - NS Equation

	strcpy(mesh->variableName, "U");
	writeResultsToFile(mesh->U,volumeVector, mesh->L_0,mesh->U_0,mesh);

	strcpy(mesh->variableName, "p");
	writeResultsToFile(mesh->p,volumeScalar, mesh->L_0,mesh->p_0,mesh);

	strcpy(mesh->variableName, "alpha");
	writeResultsToFile(mesh->alpha,volumeScalar, mesh->L_0,1.0,mesh);

	strcpy(mesh->variableName, "rho");
	writeResultsToFile(mesh->rho,volumeScalar, mesh->L_0,mesh->rho_0,mesh);

	strcpy(mesh->variableName, "mu");
	writeResultsToFile(mesh->mu,volumeScalar, mesh->L_0,mesh->mu_0,mesh);

	strcpy(mesh->variableName, "ph");
	writeResultsToFile(mesh->ph,volumeScalar, mesh->L_0,mesh->p_0,mesh);
	//=============================================================================!
	strcpy(mesh->variableName, "iMass");
	writeTimeProbe(&mesh->iMass,scalar, mesh);

	strcpy(mesh->variableName, "CPUelapsedtime");
	writeTimeProbe(&mesh->CPUelapsedtimeInSeconds,scalar, mesh);

	strcpy(mesh->variableName, "intKineticEnergy");
	writeTimeProbe(&mesh->iKE,scalar, mesh);

	strcpy(mesh->variableName, "intKineticEnergyx");
	writeTimeProbe(&mesh->iKEx,scalar, mesh);

	strcpy(mesh->variableName, "intKineticEnergyy");
	writeTimeProbe(&mesh->iKEy,scalar, mesh);

	strcpy(mesh->variableName, "intKineticEnergyz");
	writeTimeProbe(&mesh->iKEz,scalar, mesh);
	//=============================================================================!
	//	 Exact solution

	computeAndWriteAnalyticalSolution(mesh);

	// X Line Data
	//	strcpy(mesh->variableName, "hUy");
	//	writeXLineDataToFile(mesh->hUy,volumeScalar, 25, mesh->L_0,mesh->hU_0,mesh);

	// Z Line Data
	//	strcpy(mesh->variableName, "hUy");
	//	writeZLineDataToFile(mesh->hUy,volumeScalar, 50, mesh->L_0,mesh->hU_0,mesh);
	//=============================================================================!

	// X Line Data
	strcpy(mesh->variableName, "U");
	writeXLineDataToFile(mesh->U,volumeVector, 24, mesh->L_0,mesh->U_0,mesh);

	// Z Line Data
	strcpy(mesh->variableName, "U");
	writeZLineDataToFile(mesh->U,volumeVector, 24, mesh->L_0,mesh->U_0,mesh);

	//=============================================================================!
	// Constants - Maxwell Equation

	strcpy(mesh->variableName, "Sigma");
	writeResultsToFile(mesh->sigma,volumeScalar, mesh->L_0,mesh->sigma_0,mesh);

	// Primary Fields - Maxwell Equation
	strcpy(mesh->variableName, "Phi");
	writeResultsToFile(mesh->phi,volumeScalar, mesh->L_0,mesh->Phi_0,mesh);

	strcpy(mesh->variableName, "B");
	writeResultsToFile(mesh->B,volumeVector, mesh->L_0,mesh->B_0,mesh);


	// Derived Fields - Maxwell Equation

	strcpy(mesh->variableName, "Jphi");
	writeResultsToFile(mesh->Jphi,volumeVector, mesh->L_0,mesh->J_0,mesh);

	strcpy(mesh->variableName, "Ju");
	writeResultsToFile(mesh->Ju,volumeVector, mesh->L_0,mesh->J_0,mesh);

	strcpy(mesh->variableName, "J");
	writeResultsToFile(mesh->J,volumeVector, mesh->L_0,mesh->J_0,mesh);


}
