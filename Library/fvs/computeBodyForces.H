#include "mesh.h"

//#include "../Library/fvs/computeBodyForces.H"

void computeBodyForces(mesh_t *mesh);

void computeBodyForces(mesh_t *mesh){

	if (mesh->includeLorentzForce== true){

		// Lorentz Force, J cross B

		for(int cell = 0; cell < mesh->ncells; ++cell){

			dint idx	= mesh->dim_x + cell*mesh->ndimension;
			dint idy	= mesh->dim_y + cell*mesh->ndimension;
			dint idz	= mesh->dim_z + cell*mesh->ndimension;

			dfloat JcrossB_x;
			dfloat JcrossB_y;
			dfloat JcrossB_z;

			JcrossB_x = (mesh->J[idy]*mesh->B[idz] - mesh->J[idz]*mesh->B[idy]);
			JcrossB_y = (mesh->J[idz]*mesh->B[idx] - mesh->J[idx]*mesh->B[idz]);
			JcrossB_z = (mesh->J[idx]*mesh->B[idy] - mesh->J[idy]*mesh->B[idx]);

			mesh->BF3[idx]		=  mesh->alpha[cell]*JcrossB_x*mesh->Ha*mesh->Ha/mesh->Re;
			mesh->BF3[idy]		=  mesh->alpha[cell]*JcrossB_y*mesh->Ha*mesh->Ha/mesh->Re;
			mesh->BF3[idz]		=  mesh->alpha[cell]*JcrossB_z*mesh->Ha*mesh->Ha/mesh->Re;

		}


	}


	if (mesh->includeGravity== true){

		for(int cell = 0; cell < mesh->ncells; ++cell){

			dint idx	= mesh->dim_x + cell*mesh->ndimension;
			dint idy	= mesh->dim_y + cell*mesh->ndimension;
			dint idz	= mesh->dim_z + cell*mesh->ndimension;

			dfloat oneOverFrsqr = 1.0/(mesh->Fr*mesh->Fr);

			mesh->BF1[idx]		=  mesh->rho[cell]*mesh->g[idx]*oneOverFrsqr;
			mesh->BF1[idy]		=  mesh->rho[cell]*mesh->g[idy]*oneOverFrsqr;
			mesh->BF1[idz]		=  mesh->rho[cell]*mesh->g[idz]*oneOverFrsqr;

		}


		dfloat *phf;
		createField(&phf,surfaceScalar,mesh);

		interpolateCellToFaceLinear(mesh->ph,phf,volumeScalar,mesh);

		fvdGradientOfScalar(1.0, phf, mesh, mesh->BF2);

		free(phf);

	}






//	strcpy(mesh->variableName, "BF2");
//	printVolumeVector(mesh->BF2, mesh);
//
//	strcpy(mesh->variableName, "rV");
//	printVolumeScalar(mesh->rV, mesh);


}
