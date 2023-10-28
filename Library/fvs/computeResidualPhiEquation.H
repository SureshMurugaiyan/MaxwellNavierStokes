#include "mesh.h"

//#include "../Library/fvs/computeResidualPhiEquation.H"

void computeResidualPhiEquation(mesh_t *mesh);

void computeResidualPhiEquation(mesh_t *mesh){

	dfloat* sum;
	dfloat* Residual;

	createField(&sum,volumeScalar,mesh);
	createField(&Residual,volumeScalar,mesh);

	for(int cell = 0;cell < mesh->nicells; ++cell){

		sum[cell] 		= 0.0;
		Residual[cell] 	= 0.0;

		for (int face = 0; face < mesh->ncellFaces; face++) {

			int fid 	= face+cell*mesh->ncellFaces;
			int cell_n  = mesh->faceToNeighborCell[fid];

			sum[cell]  += mesh->aN[fid]*mesh->phi[cell_n];

		}

		Residual[cell] 	= mesh->bP[cell]- sum[cell] -
					  	  mesh->aP[cell]*mesh->phi[cell];


	}

	computeResidualsScalarEqn(mesh->phi,Residual,mesh);

	mesh->PhiEqnResidual_max = mesh->Residual_max ;
	mesh->PhiEqnResidual_rms = mesh->Residual_rms ;
	mesh->PhiEqnResidual_max_scaled = mesh->Residual_max_scaled ;
	mesh->PhiEqnResidual_rms_scaled = mesh->Residual_rms_scaled ;

	free(sum);
	free(Residual);


}

