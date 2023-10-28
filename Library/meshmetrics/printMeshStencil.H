#include "mesh.h"

void printMeshStencil(mesh_t *mesh);

void printMeshStencil(mesh_t *mesh){

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
	printf ( "\n" );
	printf ( "\n!--------------------------------------------- !" );
	printf ( "\n! Cell ID at different planes                  !" );
	printf ( "\n!--------------------------------------------- !" );

	for (dint k = 0; k < mesh->ncellsz; k++){

		printf("\n\n(Plane, K = %d)",k);

		if(k==0){printf("(Bottom Plane)\n");}

		if(k==mesh->ncellsz-1){printf("(Top Plane)\n");}

		printf("\n");

		for (dint i = 0; i < mesh->ncellsx; i++){

			printf("*");

			for (dint j = 0; j < mesh->ncellsy; j++){
				printf("----*");
			}

			printf("\n|");

			for (dint j = 0; j < mesh->ncellsy; j++){
				dint cell = computeLinearIndex(i,j,k, mesh->ncellsx, mesh->ncellsy,mesh->ncellsz);
				printf(" %2d |",cell);
			}

			printf("\n");
		}

		printf("*");

		for (dint j = 0; j < mesh->ncellsy; j++){
			printf("----*");
		}

	}

	printf ( "\n" );
}
