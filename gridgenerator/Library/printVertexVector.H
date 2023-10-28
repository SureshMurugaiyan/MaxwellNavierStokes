#include "mesh.h"

void printVertexVector(dfloat* Phi, mesh_t *mesh);

void printVertexVector(dfloat* Phi, mesh_t *mesh){

	for( dint k = 0; k < mesh->nverticesz; k++){

		printf("! (Plane, K = %d)\n!",k);

		for( dint i = 0; i < mesh->nverticesx; i++){
			for( dint j = 0; j < mesh->nverticesy; j++){

				dint vertex = computeLinearIndex(i,j,k, mesh->nverticesx, mesh->nverticesy,mesh->nverticesz);

				dint vid_dim_x = mesh->dim_x + vertex*mesh->ndimension ;
				dint vid_dim_y = mesh->dim_y + vertex*mesh->ndimension ;
				dint vid_dim_z = mesh->dim_z + vertex*mesh->ndimension ;


				printf(" (%.2e,%.2e,%.2e)",Phi[vid_dim_x],
										   Phi[vid_dim_y],
										   Phi[vid_dim_z]);

			}
			printf("\n!");
		}
		printf("--------------------------------------------- \n");
	}

}
