#include "mesh.h"

void meshGenerateFaceToCellConnectivity(mesh_t *mesh);

void meshCalculateInnerandBoundaryFaceID(mesh_t *mesh);
void meshFaceToOwnerCellConnectivity(mesh_t *mesh);
void meshFaceToNeighborCellConnectivity(mesh_t *mesh);

void meshCalculateInnerandBoundaryFaceIDUnitTest(mesh_t *mesh);
void meshFaceToOwnerCellConnectivityUnitTest(mesh_t *mesh);
void meshFaceToNeighborCellConnectivityUnitTest(mesh_t *mesh);

void meshGenerateFaceToCellConnectivity(mesh_t *mesh){

	meshCalculateInnerandBoundaryFaceID(mesh);

	meshFaceToOwnerCellConnectivity(mesh);

	meshFaceToNeighborCellConnectivity(mesh);

}

void meshCalculateInnerandBoundaryFaceID(mesh_t *mesh){

	dint i = 0;
	dint j = 0;

	for (dint cell = 0; cell < mesh->nicells; cell++) {

		for (dint face = 0; face < mesh->ncellFaces; face++) {

			dint fid = face + cell*mesh->ncellFaces;

			if (mesh->cellToFaces[fid] > (-1)){

				mesh->innerFaces[i] = fid;

				i++;
			}

			if (mesh->cellToFaces[fid] == (-1)){

				mesh->boundaryFaces[j] = fid;

				j++;

			}

		}

	}


//	meshCalculateInnerandBoundaryFaceIDUnitTest(mesh);

}


void meshFaceToOwnerCellConnectivity(mesh_t *mesh){

	for (dint cell = 0; cell < mesh->nicells; cell++) {
		for (dint face = 0; face < mesh->ncellFaces; face++) {

			dint fid = face+cell*mesh->ncellFaces;

			mesh->faceToOwnerCell[fid] = cell;


		}
	}


//	meshFaceToOwnerCellConnectivityUnitTest(mesh);

}

void meshFaceToNeighborCellConnectivity(mesh_t *mesh){


	for (dint cell = 0; cell < mesh->nicells; cell++) {
		for (dint face = 0; face < mesh->ncellFaces; face++) {

			dint fid 	 = face+cell*mesh->ncellFaces;
			dint cell_n  =  mesh->cellToCells[fid];

			mesh->faceToNeighborCell[fid] = cell_n;


		}
	}

	// Evaluate Neighbor cell connectivity and cell to cell connectivity for boundary cells

	int count = 0;

	for (dint cell = 0; cell < mesh->nicells; cell++) {
		for (dint face = 0; face < mesh->ncellFaces; face++) {

			dint fid 	 = face+cell*mesh->ncellFaces;
			dint cell_n  =  mesh->cellToCells[fid];

			if(cell_n == -1){

				mesh->faceToNeighborCell[fid] = mesh->nicells + count;
				mesh->cellToCells[fid] 		  = mesh->faceToNeighborCell[fid];
				count++;
			}


		}
	}


//	meshFaceToNeighborCellConnectivityUnitTest(mesh);


}


/*************************************************************/
/* UNIT TEST                                                 */
/*************************************************************/

void meshCalculateInnerandBoundaryFaceIDUnitTest(mesh_t *mesh){

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






	printf ( "\n!--------------------------------------------- !\n" );
	printf ( "! Boundary Face ID                             !\n" );
	printf ( "!--------------------------------------------- !\n" );

	for (dint face = 0; face < mesh->nBoundaryFaces; face++) {

		printf ( " %5d \n", mesh->boundaryFaces[face] );

	}

	printf ( "!--------------------------------------------- !\n" );
	printf ( "! Inner Face ID                                !\n" );
	printf ( "!--------------------------------------------- !\n" );

	for (dint face = 0; face < mesh->nInnerFaces; face++) {

		printf ( " %5d \n", mesh->innerFaces[face] );

	}
	printf ( "!--------------------------------------------- !\n" );
}

void meshFaceToOwnerCellConnectivityUnitTest(mesh_t *mesh){

	printf ( "!--------------------------------------------- !\n" );
	printf ( "! Face to OwnerCell   Connectivity             !\n" );
	printf ( "!--------------------------------------------- !\n" );

	printf ( "%11s%6s%7s%7s%6s%5s%8s\n","! Cell ID :","Left","Front","Right","Back","Top","Bottom");
	printf ( "!--------------------------------------------- !\n!" );


	for (dint cell = 0; cell < mesh->nicells; cell++) {

		printf ( " Cell%3d :", cell );

		for (dint face = 0; face < mesh->ncellFaces; face++) {

			dint fid = face 	+ cell*mesh->ncellFaces;

			printf ( "%6d", mesh->faceToOwnerCell[fid] );
		}

		printf ( "\n!" );
	}

	printf ( "--------------------------------------------- !\n" );
}


void meshFaceToNeighborCellConnectivityUnitTest(mesh_t *mesh){

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

	printf ( "\n!--------------------------------------------- !\n" );
	printf ( "! Face to NeighborCell Connectivity            !\n" );
	printf ( "!--------------------------------------------- !\n" );

	printf ( "%11s%6s%7s%7s%6s%5s%8s\n","! Cell ID :","Left","Front","Right","Back","Top","Bottom");
	printf ( "!--------------------------------------------- !\n!" );


	for (dint cell = 0; cell < mesh->nicells; cell++) {

		printf ( " Cell%3d :", cell );

		for (dint face = 0; face < mesh->ncellFaces; face++) {

			dint fid = face 	+ cell*mesh->ncellFaces;

			printf ( "%6d", mesh->faceToNeighborCell[fid] );
		}

		printf ( "\n!" );
	}

	printf ( "--------------------------------------------- !\n" );

	printf ( "!--------------------------------------------- !\n" );
	printf ( "! Cell To Cell Connectivity                    !\n" );
	printf ( "!--------------------------------------------- !\n" );

	printf ( "%11s%6s%7s%7s%6s%5s%8s\n","! Cell ID :","Left","Front","Right","Back","Top","Bottom");
	printf ( "!--------------------------------------------- !\n!" );

	for (dint cell = 0; cell < mesh->nicells; cell++) {

		printf ( " Cell%3d :", cell );

		for (dint face = 0; face < mesh->ncellFaces; face++) {

			dint fid = face 	+ cell*mesh->ncellFaces;

			printf ( "%6d", mesh->cellToCells[fid] );
		}

		printf ( "\n!" );
	}

	printf ( "--------------------------------------------- !\n" );
}
