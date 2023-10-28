#include "mesh.h"

//#include "../Library/fvs/computeResidualsVectorEqn.H"

void computeResidualsVectorEqn(dfloat* phi,dfloat* Residual,mesh_t *mesh);

void computeResidualsVectorEqn(dfloat* phi,dfloat* Residual,mesh_t *mesh){


	dint cell;
	dint dimx; // dim = 0 gives x component
	dint dimy; // dim = 1 gives y component
	dint dimz; // dim = 2 gives z component
	dint idx;
	dint idy;
	dint idz;



	dfloat* magResidual;
	createField(&magResidual,volumeVector,mesh);


	// 1. Maximum Residual

	for (int cell = 0; cell < mesh->nicells; cell++) {
		for (int dim = 0; dim < mesh->ndimension; dim++) {

			dint cell_dim = dim + cell*mesh->ndimension;

			magResidual[cell_dim] = magnitude(Residual[cell_dim]);

		}
	}

	// initialize

	cell = 0;

	dimx = 0; // dim = 0 gives x component
	dimy = 1; // dim = 1 gives y component
	dimz = 2; // dim = 2 gives z component

	idx	= dimx+cell*mesh->ndimension;
	idy	= dimy+cell*mesh->ndimension;
	idz	= dimz+cell*mesh->ndimension;

	mesh->Residual_max_Xeqn = magResidual[idx];
	mesh->Residual_max_Yeqn = magResidual[idy];
	mesh->Residual_max_Zeqn = magResidual[idz];

	for (int cell = 0; cell < mesh->nicells; cell++) {

		dimx = 0; // dim = 0 gives x component
		dimy = 1; // dim = 1 gives y component
		dimz = 2; // dim = 2 gives z component

		idx	= dimx+cell*mesh->ndimension;
		idy	= dimy+cell*mesh->ndimension;
		idz	= dimz+cell*mesh->ndimension;

		mesh->Residual_max_Xeqn = mymax(magResidual[idx],mesh->Residual_max_Xeqn);
		mesh->Residual_max_Yeqn = mymax(magResidual[idy],mesh->Residual_max_Yeqn);
		mesh->Residual_max_Zeqn = mymax(magResidual[idz],mesh->Residual_max_Zeqn);
	}

	// 2. RMS Residual

	dfloat sum_x = 0;
	dfloat sum_y = 0;
	dfloat sum_z = 0;

	for (int cell = 0; cell < mesh->nicells; cell++) {

		dimx = 0; // dim = 0 gives x component
		dimy = 1; // dim = 1 gives y component
		dimz = 2; // dim = 2 gives z component

		idx	= dimx+cell*mesh->ndimension;
		idy	= dimy+cell*mesh->ndimension;
		idz	= dimz+cell*mesh->ndimension;


		sum_x = sum_x + magResidual[idx]*magResidual[idx];
		sum_y = sum_y + magResidual[idy]*magResidual[idy];
		sum_z = sum_z + magResidual[idz]*magResidual[idz];
	}

	mesh->Residual_rms_Xeqn = sqrt(sum_x/mesh->nicells) ;
	mesh->Residual_rms_Yeqn = sqrt(sum_y/mesh->nicells) ;
	mesh->Residual_rms_Zeqn = sqrt(sum_z/mesh->nicells) ;

	// Scaled Residuals

	dfloat scalefactor_x;
	dfloat scalefactor_y;
	dfloat scalefactor_z;

	dfloat phi_scale_x;
	dfloat phi_scale_y;
	dfloat phi_scale_z;

	dfloat phi_max_x;
	dfloat phi_max_y;
	dfloat phi_max_z;

	dfloat phi_min_x;
	dfloat phi_min_y;
	dfloat phi_min_z;

	dfloat del_phi_max_x;
	dfloat del_phi_max_y;
	dfloat del_phi_max_z;

	dfloat abs_phi_max_x;
	dfloat abs_phi_max_y;
	dfloat abs_phi_max_z;

	// initialize

	cell = 0;

	dimx = 0; // dim = 0 gives x component
	dimy = 1; // dim = 1 gives y component
	dimz = 2; // dim = 2 gives z component

	idx	= dimx+cell*mesh->ndimension;
	idy	= dimy+cell*mesh->ndimension;
	idz	= dimz+cell*mesh->ndimension;

	phi_max_x = phi[idx];
	phi_max_y = phi[idy];
	phi_max_z = phi[idz];

	phi_min_x = phi[idx];
	phi_min_y = phi[idy];
	phi_min_z = phi[idz];

	for (int cell = 0; cell < mesh->nicells; cell++) {

		dimx = 0; // dim = 0 gives x component
		dimy = 1; // dim = 1 gives y component
		dimz = 2; // dim = 2 gives z component

		idx	= dimx+cell*mesh->ndimension;
		idy	= dimy+cell*mesh->ndimension;
		idz	= dimz+cell*mesh->ndimension;

		phi_max_x = mymax(phi[idx],phi_max_x);
		phi_max_y = mymax(phi[idy],phi_max_y);
		phi_max_z = mymax(phi[idz],phi_max_z);

		phi_min_x = mymin(phi[idx],phi_min_x);
		phi_min_y = mymin(phi[idy],phi_min_y);
		phi_min_z = mymin(phi[idz],phi_min_z);


	}

	del_phi_max_x = phi_max_x - phi_min_x;
	del_phi_max_y = phi_max_y - phi_min_y;
	del_phi_max_z = phi_max_z - phi_min_z;

	abs_phi_max_x = magnitude(phi_max_x);
	abs_phi_max_y = magnitude(phi_max_y);
	abs_phi_max_z = magnitude(phi_max_z);

	phi_scale_x = mymax(del_phi_max_x,abs_phi_max_x);
	phi_scale_y = mymax(del_phi_max_y,abs_phi_max_y);
	phi_scale_z = mymax(del_phi_max_z,abs_phi_max_z);


	for (int cell = 0; cell < mesh->nicells; cell++) {

		dimx = 0; // dim = 0 gives x component
		dimy = 1; // dim = 1 gives y component
		dimz = 2; // dim = 2 gives z component

		idx	= dimx+cell*mesh->ndimension;
		idy	= dimy+cell*mesh->ndimension;
		idz	= dimz+cell*mesh->ndimension;

		scalefactor_x = mesh->AP[idx]*phi_scale_x;
		scalefactor_y = mesh->AP[idy]*phi_scale_y;
		scalefactor_z = mesh->AP[idz]*phi_scale_z;

		magResidual[idx] = magResidual[idx] /stabilise(scalefactor_x,mesh);
		magResidual[idy] = magResidual[idy] /stabilise(scalefactor_y,mesh);
		magResidual[idz] = magResidual[idz] /stabilise(scalefactor_z,mesh);
	}


	// 3. Maximum Residual Scaled


	// initialize

	cell = 0;

	dimx = 0; // dim = 0 gives x component
	dimy = 1; // dim = 1 gives y component
	dimz = 2; // dim = 2 gives z component

	idx	= dimx+cell*mesh->ndimension;
	idy	= dimy+cell*mesh->ndimension;
	idz	= dimz+cell*mesh->ndimension;

	mesh->Residual_max_scaled_Xeqn = magResidual[idx];
	mesh->Residual_max_scaled_Yeqn = magResidual[idy];
	mesh->Residual_max_scaled_Zeqn = magResidual[idz];


	for (int cell = 0; cell < mesh->nicells; cell++) {

		dimx = 0; // dim = 0 gives x component
		dimy = 1; // dim = 1 gives y component
		dimz = 2; // dim = 2 gives z component

		idx	= dimx+cell*mesh->ndimension;
		idy	= dimy+cell*mesh->ndimension;
		idz	= dimz+cell*mesh->ndimension;

		mesh->Residual_max_scaled_Xeqn = mymax(magResidual[idx],mesh->Residual_max_scaled_Xeqn);
		mesh->Residual_max_scaled_Yeqn = mymax(magResidual[idy],mesh->Residual_max_scaled_Yeqn);
		mesh->Residual_max_scaled_Zeqn = mymax(magResidual[idz],mesh->Residual_max_scaled_Zeqn);
	}


	// 4. RMS Residual Scaled


	sum_x = 0;
	sum_y = 0;
	sum_z = 0;


	for (int cell = 0; cell < mesh->nicells; cell++) {

		dimx = 0; // dim = 0 gives x component
		dimy = 1; // dim = 1 gives y component
		dimz = 2; // dim = 2 gives z component

		idx	= dimx+cell*mesh->ndimension;
		idy	= dimy+cell*mesh->ndimension;
		idz	= dimz+cell*mesh->ndimension;

		sum_x = sum_x + magResidual[idx]*magResidual[idx];
		sum_y = sum_y + magResidual[idy]*magResidual[idy];
		sum_z = sum_z + magResidual[idz]*magResidual[idz];
	}

	mesh->Residual_rms_scaled_Xeqn = sqrt(sum_x/mesh->nicells) ;
	mesh->Residual_rms_scaled_Yeqn = sqrt(sum_y/mesh->nicells) ;
	mesh->Residual_rms_scaled_Zeqn = sqrt(sum_z/mesh->nicells) ;


	free(magResidual);

}
