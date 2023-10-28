#include "mesh.h"
#include "../../setUpCase/FieldNames.H"

//#include "../Library/fvs/constructUEquation.H"

typedef struct{

	dfloat *gradP;     // Pressure gradient


}ueqn_t;



void constructUEquation(mesh_t *mesh);

void computeDerivedFieldsUEquation(ueqn_t *ueqn,mesh_t *mesh);

void constructUEquation(mesh_t *mesh){

	ueqn_t *ueqn = (ueqn_t*) calloc(1, sizeof(ueqn_t));

	createField(&ueqn->gradP,volumeVector,mesh);


	computeDerivedFieldsUEquation(ueqn,mesh);

	reinitializeField(mesh->AP,volumeVector,mesh);
	reinitializeField(mesh->AN,surfaceVector,mesh);
	reinitializeField(mesh->BP,volumeVector,mesh);


	dfloat Ueqn	=	0.0;


	// Time derivative
	Ueqn = fvmddtVectorWithScalarCoefficient(1.0,mesh->rho,mesh->U, U, euler,mesh);

	// Laplacian of U

	dfloat oneOverRe = 1.0/mesh->Re;

	Ueqn = fvmLaplacianOfVectorWithScalarCoefficient(-oneOverRe,mesh->muf, mesh->U,U,mesh);





	// Divergence of U U
	mesh->divScheme				= Gamma;
	mesh->divSchemeblendFactor	= 0.2;

	Ueqn = fvmDivergenceOfVector(1.0, mesh->U,U,mesh->rhophiU,mesh->divScheme,mesh->divSchemeblendFactor,mesh);

	if (mesh->includeGravity== true){

		Ueqn = fvcVolumetricSourceVector(-1.0,mesh->BF1,mesh);

		Ueqn = fvcVolumetricSourceVector(1.0,mesh->BF2,mesh);

	}

	if (mesh->includeLorentzForce== true){

		Ueqn = fvcVolumetricSourceVector(-1.0,mesh->BF3,mesh);

	}





	storeUEqnMatrixCoefficients(mesh);

	Ueqn = fvcVolumetricSourceVector(1.0,ueqn->gradP,mesh);





	free(ueqn->gradP);

	free(ueqn);

}

void computeDerivedFieldsUEquation(ueqn_t *ueqn,mesh_t *mesh){

	dfloat *pf;
	createField(&pf,surfaceScalar,mesh);

	interpolateCellToFaceLinear(mesh->p,pf,volumeScalar,mesh);

	fvdGradientOfScalar(1.0, pf, mesh, ueqn->gradP);

	free(pf);

}
