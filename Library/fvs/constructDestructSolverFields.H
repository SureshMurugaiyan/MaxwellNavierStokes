#include "mesh.h"

//#include "../Library/fvs/constructDestructSolverFields.H"

void constructSolverFields(mesh_t *mesh);
void destructSolverFields(mesh_t *mesh);

void constructSolverFields(mesh_t *mesh){

	// Calculating the Matrix coefficients
	createField(&mesh->aP,volumeScalar,mesh);
	createField(&mesh->aN,surfaceScalar,mesh);
	createField(&mesh->bP,volumeScalar,mesh);

	createField(&mesh->AP,volumeVector,mesh);
	createField(&mesh->AN,surfaceVector,mesh);
	createField(&mesh->BP,volumeVector,mesh);

	// For Computing HbyA, store coefficients without grad P term

	createField(&mesh->AAP,volumeVector,mesh);
	createField(&mesh->AAN,surfaceVector,mesh);
	createField(&mesh->BBP,volumeVector,mesh);
}

void destructSolverFields(mesh_t *mesh){

	free(mesh->aP);
	free(mesh->aN);
	free(mesh->bP);

	free(mesh->AP);
	free(mesh->AN);
	free(mesh->BP);

	free(mesh->AAP);
	free(mesh->AAN);
	free(mesh->BBP);


}

