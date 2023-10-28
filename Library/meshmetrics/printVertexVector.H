
#include "mesh.h"

void printVertexVector(dfloat* Phi, mesh_t *mesh);

void printVertexVector(dfloat* Phi, mesh_t *mesh){


	printf ( "\n!--------------------------------------------- !" );
	printf ( "\n! Geometry                                     !" );
	printf ( "\n!--------------------------------------------- !" );
	printf ( "\n*"
	"\n*                ^ Z(k)"
	"\n*                ."
	"\n*                ."
    "\n*                ."
    "\n*               B*--------------------------*D"
	"\n*               *                          *|"
    "\n*              * .                        * |"
	"\n*             *  .                       *  |"
	"\n*            *   .                      *   |"
	"\n*           *    .                     *    |"
	"\n*          *     .                    *     |"
	"\n*        F*--------------------------*H     |"
	"\n*         |      .                   |      |"
	"\n*         |      .                   |      |"
	"\n*         |     A. . . . . . . . . . | . . .*C-----------> Y(j)"
	"\n*         |     .                    |     *"
	"\n*         |   .                      |   *"
	"\n*         | .                        | *"
	"\n*        E*--------------------------*G"
	"\n*       ."
	"\n*     ."
	"\n*   ."
	"\n*  X(i)");

	printf ( "\n!--------------------------------------------- !" );
	printf ( "\n!%3s                                           !",mesh->variableName);
	printf ( "\n!--------------------------------------------- !" );


	for( dint k = 0; k < mesh->ncellsz; k++){

		printf("\n! (Plane, K = %d)",k);

		for( dint i = 0; i < mesh->ncellsx; i++){

			printf("\n!");

			for( dint j = 0; j < mesh->ncellsy; j++){

				dint cell = computeLinearIndex(i,j,k, mesh->ncellsx, mesh->ncellsy,mesh->ncellsz);

				dint vid_vertA = mesh->vertA + cell*mesh->ncellVertices;
				dint vid_vertC = mesh->vertC + cell*mesh->ncellVertices;

				dint vid_vertA_x = mesh->dim_x + vid_vertA*mesh->ndimension;
				dint vid_vertA_y = mesh->dim_y + vid_vertA*mesh->ndimension;
				dint vid_vertA_z = mesh->dim_z + vid_vertA*mesh->ndimension;

				dint vid_vertC_x = mesh->dim_x + vid_vertC*mesh->ndimension;
				dint vid_vertC_y = mesh->dim_y + vid_vertC*mesh->ndimension;
				dint vid_vertC_z = mesh->dim_z + vid_vertC*mesh->ndimension;

				// Print A vertices
				printf ( " (%9.2e,%9.2e,%9.2e)",Phi[vid_vertA_x],Phi[vid_vertA_y],Phi[vid_vertA_z]);

				// Print C vertices for last column
				if (j==(mesh->ncellsy-1)){
					printf ( " (%9.2e,%9.2e,%9.2e)",Phi[vid_vertC_x],Phi[vid_vertC_y],Phi[vid_vertC_z] );
				}
			}
		}


		printf("\n!");

		for( dint j = 0; j < mesh->ncellsy; j++){

			dint i = (mesh->ncellsx-1);

			dint cell = computeLinearIndex(i,j,k, mesh->ncellsx, mesh->ncellsy,mesh->ncellsz);

			dint vid_vertE = mesh->vertE + cell*mesh->ncellVertices;
			dint vid_vertG = mesh->vertG + cell*mesh->ncellVertices;

			dint vid_vertE_x = mesh->dim_x + vid_vertE*mesh->ndimension;
			dint vid_vertE_y = mesh->dim_y + vid_vertE*mesh->ndimension;
			dint vid_vertE_z = mesh->dim_z + vid_vertE*mesh->ndimension;

			dint vid_vertG_x = mesh->dim_x + vid_vertG*mesh->ndimension;
			dint vid_vertG_y = mesh->dim_y + vid_vertG*mesh->ndimension;
			dint vid_vertG_z = mesh->dim_z + vid_vertG*mesh->ndimension;

			// Print E vertices for last row
			printf ( " (%9.2e,%9.2e,%9.2e)",Phi[vid_vertE_x],Phi[vid_vertE_y],Phi[vid_vertE_z]);

			// Print G vertices for last row and last column
			if  (j==(mesh->ncellsy-1)){
				printf ( " (%9.2e,%9.2e,%9.2e)",Phi[vid_vertG_x],Phi[vid_vertG_y],Phi[vid_vertG_z] );
			}

		}

		printf("\n!---------------------------------------------!");

	}


	// Print for last top plane

	dint k = (mesh->ncellsz-1);

	printf("\n! (Plane, K = %d)",k+1);

	for( dint i = 0; i < mesh->ncellsx; i++){

		printf("\n!");

		for( dint j = 0; j < mesh->ncellsy; j++){

			dint cell = computeLinearIndex(i,j,k, mesh->ncellsx, mesh->ncellsy,mesh->ncellsz);

			dint vid_vertB = mesh->vertB + cell*mesh->ncellVertices;
			dint vid_vertD = mesh->vertD + cell*mesh->ncellVertices;

			dint vid_vertB_x = mesh->dim_x + vid_vertB*mesh->ndimension;
			dint vid_vertB_y = mesh->dim_y + vid_vertB*mesh->ndimension;
			dint vid_vertB_z = mesh->dim_z + vid_vertB*mesh->ndimension;

			dint vid_vertD_x = mesh->dim_x + vid_vertD*mesh->ndimension;
			dint vid_vertD_y = mesh->dim_y + vid_vertD*mesh->ndimension;
			dint vid_vertD_z = mesh->dim_z + vid_vertD*mesh->ndimension;

			// Print B vertices
			printf ( " (%9.2e,%9.2e,%9.2e)",Phi[vid_vertB_x],Phi[vid_vertB_y],Phi[vid_vertB_z] );

			// Print D vertices for last column
			if (j==(mesh->ncellsy-1)){
				printf ( " (%9.2e,%9.2e,%9.2e)",Phi[vid_vertD_x],Phi[vid_vertD_y],Phi[vid_vertD_z] );
			}

		}
	}

	printf("\n!");

	for( dint j = 0; j < mesh->ncellsy; j++){

		dint i = (mesh->ncellsx-1);

		dint cell = computeLinearIndex(i,j,k, mesh->ncellsx, mesh->ncellsy,mesh->ncellsz);

		dint vid_vertF = mesh->vertF + cell*mesh->ncellVertices;
		dint vid_vertH = mesh->vertH + cell*mesh->ncellVertices;

		dint vid_vertF_x = mesh->dim_x + vid_vertF*mesh->ndimension;
		dint vid_vertF_y = mesh->dim_y + vid_vertF*mesh->ndimension;
		dint vid_vertF_z = mesh->dim_z + vid_vertF*mesh->ndimension;

		dint vid_vertH_x = mesh->dim_x + vid_vertH*mesh->ndimension;
		dint vid_vertH_y = mesh->dim_y + vid_vertH*mesh->ndimension;
		dint vid_vertH_z = mesh->dim_z + vid_vertH*mesh->ndimension;

		// Print south west vertices for last row
		printf ( " (%9.2e,%9.2e,%9.2e)",Phi[vid_vertF_x],Phi[vid_vertF_y],Phi[vid_vertF_z]);

		// Print south east vertices for last row last column
		if  (j==(mesh->ncellsy-1)){
			printf ( " (%9.2e,%9.2e,%9.2e)",Phi[vid_vertH_x],Phi[vid_vertH_y],Phi[vid_vertH_z]);
		}

	}

	printf("\n!---------------------------------------------!");

}
