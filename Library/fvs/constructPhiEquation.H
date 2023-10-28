#include "mesh.h"

#include "../../setUpCase/FieldNames.H"

//#include "../Library/fvs/constructPhiEquation.H"


void constructPhiEquation(mesh_t *mesh);

void computeDerivedFieldsPhiEquation(mesh_t *mesh);

void constructPhiEquation(mesh_t *mesh){


	computeDerivedFieldsPhiEquation(mesh);

	reinitializeField(mesh->aP,volumeScalar,mesh);
	reinitializeField(mesh->aN,surfaceScalar,mesh);
	reinitializeField(mesh->bP,volumeScalar,mesh);

	dint Phieqn;

	Phieqn = fvmLaplacianOfScalarWithScalarCoefficient(1.0, mesh->sigmaf,mesh->phi,phi,mesh);

//	Phieqn = fvcDivergenceOfFlux (-1.0, mesh->sigmaUcrossBFlux,mesh);

}

void computeDerivedFieldsPhiEquation(mesh_t *mesh){

	dfloat *sigmaUcrossBf;
	createField(&sigmaUcrossBf,surfaceVector,mesh);


	interpolateCellToFaceLinear(mesh->sigma,mesh->sigmaf,volumeScalar,mesh);


	for (int cell = 0; cell < mesh->ncells; cell++) {

		int idx,idy,idz;
		dfloat UXB_x,UXB_y,UXB_z;

		idx	= mesh->dim_x+cell*mesh->ndimension;
		idy	= mesh->dim_y+cell*mesh->ndimension;
		idz	= mesh->dim_z+cell*mesh->ndimension;

		UXB_x =   (mesh->U[idy]*mesh->B[idz] - mesh->U[idz]*mesh->B[idy]);
		UXB_y =   (mesh->U[idz]*mesh->B[idx] - mesh->U[idx]*mesh->B[idz]);
		UXB_z =   (mesh->U[idx]*mesh->B[idy] - mesh->U[idy]*mesh->B[idx]);

		mesh->Ju[idx] = mesh->sigma[cell]*UXB_x;
		mesh->Ju[idy] = mesh->sigma[cell]*UXB_y;
		mesh->Ju[idz] = mesh->sigma[cell]*UXB_z;

	}

	interpolateCellToFaceLinear(mesh->Ju,sigmaUcrossBf,volumeVector,mesh);

	constructFluxofSurfaceVector(sigmaUcrossBf,mesh->sigmaUcrossBFlux,mesh);



	//	strcpy(mesh->variableName, "sigmaf");
	//	printSurfaceScalar(mesh->sigmaf, mesh);


	free(sigmaUcrossBf);


}
