#include "mesh.h"

//#include "../Library/fvs/computeResidualsScalarEqn.H"


void computeResidualsScalarEqn(dfloat* phi,dfloat* Residual,mesh_t *mesh);

void computeResidualsScalarEqn(dfloat* phi,dfloat* Residual,mesh_t *mesh){

	dfloat* magResidual;
	createField(&magResidual,volumeScalar,mesh);


	// 1. Maximum Residual

	for (int cell = 0; cell < mesh->nicells; cell++) {

		magResidual[cell] = magnitude(Residual[cell]);
	}

	mesh->Residual_max = magResidual[0];

	for (int cell = 0; cell < mesh->nicells; cell++) {

		mesh->Residual_max = mymax(magResidual[cell],mesh->Residual_max);
	}

	// 2. RMS Residual

	dfloat sum = 0;

	for (int cell = 0; cell < mesh->nicells; cell++) {

		sum = sum + magResidual[cell]*magResidual[cell];
	}


	mesh->Residual_rms = sqrt(sum/mesh->nicells) ;

	// Scaled Residuals

	dfloat scalefactor;
	dfloat phi_scale;
	dfloat phi_max;
	dfloat phi_min;
	dfloat del_phi_max;
	dfloat abs_phi_max;

	phi_max = phi[0];
	phi_min = phi[0];

	for (int cell = 0; cell < mesh->nicells; cell++) {

		phi_max = mymax(phi[cell],phi_max);
		phi_min = mymin(phi[cell],phi_min);
	}

	del_phi_max = phi_max - phi_min;
	abs_phi_max = magnitude(phi_max);


	phi_scale = mymax(del_phi_max,abs_phi_max);


	for (int cell = 0; cell < mesh->nicells; cell++) {

		scalefactor = mesh->aP[cell]*phi_scale;

		magResidual[cell] = magResidual[cell] /stabilise(scalefactor,mesh);
	}

	// 3. Maximum Residual Scaled

	mesh->Residual_max_scaled = magResidual[0];

	for (int cell = 0; cell < mesh->nicells; cell++) {

		mesh->Residual_max_scaled = mymax(magResidual[cell],mesh->Residual_max_scaled);
	}

	// 4. RMS Residual Scaled

	sum = 0;

	for (int cell = 0; cell < mesh->nicells; cell++) {

		sum = sum + magResidual[cell]*magResidual[cell];
	}

	mesh->Residual_rms_scaled = sqrt(sum/mesh->nicells) ;



	free(magResidual);
}
