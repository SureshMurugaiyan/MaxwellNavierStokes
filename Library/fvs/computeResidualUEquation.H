#include "mesh.h"

//#include "../Library/fvs/computeResidualUEquation.H"


void computeResidualUEquation(mesh_t *mesh);

void computeResidualUEquation(mesh_t *mesh){

	dfloat* sum;
	dfloat* Residual;

	createField(&sum,volumeVector,mesh);
	createField(&Residual,volumeVector,mesh);

	for(int cell = 0;cell < mesh->nicells; ++cell){

		for (int dim = 0; dim < mesh->ndimension; dim++) {

			dint cell_dim = dim + cell*mesh->ndimension;

			sum[cell_dim] 		= 0.0;
			Residual[cell_dim] 	= 0.0;
		}

		for (int face = 0; face < mesh->ncellFaces; face++) {

			dint fid 	= face+cell*mesh->ncellFaces;
			dint cell_n	=  mesh->faceToNeighborCell[fid];

			for (dint dim = 0; dim < mesh->ndimension; dim++) {

				dint cell_dim 	= dim + cell*mesh->ndimension;
				dint fid_dim  	= dim + fid*mesh->ndimension;
				dint cell_ndim 	= dim + cell_n*mesh->ndimension;

				sum[cell_dim] += (mesh->AN[fid_dim]*mesh->U[cell_ndim]);

			}
		}

		for (int dim = 0; dim < mesh->ndimension; dim++) {

			dint cell_dim = dim + cell*mesh->ndimension;


			Residual[cell_dim] = (mesh->BP[cell_dim] - sum[cell_dim] -
								  mesh->AP[cell_dim]*mesh->U[cell_dim]);
		}

	}

	computeResidualsVectorEqn(mesh->U,Residual,mesh);

	mesh->UEqnResidual_max_Xeqn = mesh->Residual_max_Xeqn ;
	mesh->UEqnResidual_max_Yeqn = mesh->Residual_max_Yeqn ;
	mesh->UEqnResidual_max_Zeqn = mesh->Residual_max_Zeqn ;

	mesh->UEqnResidual_rms_Xeqn = mesh->Residual_rms_Xeqn ;
	mesh->UEqnResidual_rms_Yeqn = mesh->Residual_rms_Yeqn ;
	mesh->UEqnResidual_rms_Zeqn = mesh->Residual_rms_Zeqn ;

	mesh->UEqnResidual_max_scaled_Xeqn = mesh->Residual_max_scaled_Xeqn ;
	mesh->UEqnResidual_max_scaled_Yeqn = mesh->Residual_max_scaled_Yeqn ;
	mesh->UEqnResidual_max_scaled_Zeqn = mesh->Residual_max_scaled_Zeqn ;

	mesh->UEqnResidual_rms_scaled_Xeqn = mesh->Residual_rms_scaled_Xeqn ;
	mesh->UEqnResidual_rms_scaled_Yeqn = mesh->Residual_rms_scaled_Yeqn ;
	mesh->UEqnResidual_rms_scaled_Zeqn = mesh->Residual_rms_scaled_Zeqn ;


	free(sum);
	free(Residual);


}

