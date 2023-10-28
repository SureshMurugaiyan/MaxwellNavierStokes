#include "mesh.h"

void meshGenerateCellConnectivity(mesh_t *mesh);

void meshCellToVerticesConnectivity(mesh_t *mesh);
void meshCellToCellsConnectivity(mesh_t *mesh);
void meshCellToFacesConnectivity(mesh_t *mesh);
void meshFaceToVerticesConnectivity(mesh_t *mesh);

void meshCellToVerticesConnectivityUnitTest(mesh_t *mesh);
void meshCellToCellsConnectivityUnitTest(mesh_t *mesh);
void meshCellToFacesConnectivityUnitTest(mesh_t *mesh);
void meshFaceToVerticesConnectivityUnitTest(mesh_t *mesh);


void meshGenerateCellConnectivity(mesh_t *mesh){

	meshCellToVerticesConnectivity(mesh);

	meshCellToCellsConnectivity(mesh);

	meshCellToFacesConnectivity(mesh);

	meshFaceToVerticesConnectivity(mesh);


}

void meshCellToVerticesConnectivity(mesh_t *mesh){


	for(dint cell = 0; cell < mesh->nicells; cell++){
		for(dint vert = 0; vert < mesh->ncellVertices; vert++){

			dint vid = vert + cell*mesh->ncellVertices;

			mesh->cellToVertex[vid] = vid;

		}
	}

//	meshCellToVerticesConnectivityUnitTest(mesh);

}



void meshCellToCellsConnectivity(mesh_t *mesh){


	/*                  ^ Z(k)
	 *                  .
	 *                  .
     *                  .
     *                 B*--------------------------*D
	 *				   *.                         *|
     *                * .                        * |
	 *      	     *  .                       *  |
	 *  	  	    *   .                      *   |
	 *  	  	   *    .                     *    |
	 * 	          *		.	                 *     |
	 *          F*--------------------------*H     |
	 *  	     |      .                   |      |
	 * 	     	 |      .                   |      |
	 *           |	   A. . . . . . . . . . | . . .*C-----------> Y(j)
	 *           |     .                    |     *
	 *       	 |   .                      |   *
	 *  	     | .                        | *
	 *          E*--------------------------*G
	 * 		  .
	        .
	      .
	     X(i)
	*/


	dint *cells;
	cells	= (dint*) calloc(mesh->ncellFaces, sizeof(dint));

	for( dint i = 0; i < mesh->ncellsx; i++){
		for( dint j = 0; j < mesh->ncellsy; j++){
			for( dint k = 0; k < mesh->ncellsz; k++){

				dint cell = computeLinearIndex(i,j,k, mesh->ncellsx,mesh->ncellsy,mesh->ncellsz);

				cells[0] =  computeLinearIndex(i,j-1,k, mesh->ncellsx,mesh->ncellsy,mesh->ncellsz);  // faceAEFB Neighbor
				cells[1] =  computeLinearIndex(i+1,j,k, mesh->ncellsx,mesh->ncellsy,mesh->ncellsz);  // faceEGHF Neighbor
				cells[2] =  computeLinearIndex(i,j+1,k, mesh->ncellsx,mesh->ncellsy,mesh->ncellsz);  // faceGCDH Neighbor
				cells[3] =  computeLinearIndex(i-1,j,k, mesh->ncellsx,mesh->ncellsy,mesh->ncellsz);  // faceCABD Neighbor
				cells[4] =  computeLinearIndex(i,j,k+1, mesh->ncellsx,mesh->ncellsy,mesh->ncellsz);  // faceFHDB Neighbor
				cells[5] =  computeLinearIndex(i,j,k-1, mesh->ncellsx,mesh->ncellsy,mesh->ncellsz);  // faceACGE Neighbor


				for (dint face = 0; face < mesh->ncellFaces; face++) {

					dint fid = face + cell*mesh->ncellFaces;

					mesh->cellToCells[fid] = cells[face];


					// Corrections for Boundary

					// Correct for faceAEFB boundary
					if ((j==0) && (face==mesh->faceAEFB)){
						mesh->cellToCells[fid] = -1;
					}

					// Correct for faceGCDH boundary
					if ((j ==(mesh->ncellsy-1)) && (face==mesh->faceGCDH)){
						mesh->cellToCells[fid] = -1;
					}


					// Correct for faceEGHF boundary
					if ((i==(mesh->ncellsx-1)) && (face==mesh->faceEGHF)){
						mesh->cellToCells[fid] = -1;
					}

					// Correct for faceCABD boundary
					if ((i==0) && (face==mesh->faceCABD)){
						mesh->cellToCells[fid] = -1;
					}

					// Correct for faceFHDB boundary
					if ((k==(mesh->ncellsz-1)) && (face==mesh->faceFHDB)){
						mesh->cellToCells[fid] = -1;
					}

					// Correct for faceACGE boundary
					if ((k==0) && (face==mesh->faceACGE)){
						mesh->cellToCells[fid] = -1;
					}

				}





			}
		}
	}

	free(cells);


//	meshCellToCellsConnectivityUnitTest(mesh);



}


void meshCellToFacesConnectivity(mesh_t *mesh){

	/*                  ^ Z(k)
	 *                  .
	 *                  .
     *                  .
     *                 B*--------------------------*D
	 *				   *.                         *|
     *                * .                        * |
	 *      	     *  .                       *  |
	 *  	  	    *   .                      *   |
	 *  	  	   *    .                     *    |
	 * 	          *		.	                 *     |
	 *          F*--------------------------*H     |
	 *  	     |      .                   |      |
	 * 	     	 |      .                   |      |
	 *           |	   A. . . . . . . . . . | . . .*C-----------> Y(j)
	 *           |     .                    |     *
	 *       	 |   .                      |   *
	 *  	     | .                        | *
	 *          E*--------------------------*G
	 * 		  .
	        .
	      .
	     X(i)
	*/


	dint *neighborfaces;
	neighborfaces	= (dint*) calloc(mesh->ncellFaces, sizeof(dint));

	neighborfaces[mesh->faceAEFB] =  mesh->faceGCDH;
	neighborfaces[mesh->faceEGHF] =  mesh->faceCABD;
	neighborfaces[mesh->faceGCDH] =  mesh->faceAEFB;
	neighborfaces[mesh->faceCABD] =  mesh->faceEGHF;
	neighborfaces[mesh->faceFHDB] =  mesh->faceACGE;
	neighborfaces[mesh->faceACGE] =  mesh->faceFHDB;


	for(dint cell = 0; cell < mesh->nicells; cell++){

		for (dint face = 0; face < mesh->ncellFaces; face++) {

			dint fid = face + cell*mesh->ncellFaces;

			mesh->cellToFaces[fid] = neighborfaces[face];

		}

	}

	free(neighborfaces);



	for (dint cell = 0; cell < mesh->nicells; cell++) {

		for (dint face = 0; face < mesh->ncellFaces; face++) {

			dint fid    = face + cell*mesh->ncellFaces;
			dint cell_n = mesh->cellToCells[fid];

			if ( cell_n==(-1)){

				mesh->cellToFaces[fid] = -1;
			}

		}

	}





//	meshCellToFacesConnectivityUnitTest(mesh);

}

void meshFaceToVerticesConnectivity(mesh_t *mesh){

	/*                  ^ Z(k)
	 *                  .
	 *                  .
     *                  .
     *                 B*--------------------------*D
	 *				   *.                         *|
     *                * .                        * |
	 *      	     *  .                       *  |
	 *  	  	    *   .                      *   |
	 *  	  	   *    .                     *    |
	 * 	          *		.	                 *     |
	 *          F*--------------------------*H     |
	 *  	     |      .                   |      |
	 * 	     	 |      .                   |      |
	 *           |	   A. . . . . . . . . . | . . .*C-----------> Y(j)
	 *           |     .                    |     *
	 *       	 |   .                      |   *
	 *  	     | .                        | *
	 *          E*--------------------------*G
	 * 		  .
	        .
	      .
	     X(i)
	*/

	/*                  ^ Z(k)
	 *                  .
	 *                  .
     *                  .
     *                1 *--------------------------*
	 *				   *.                         *|3
     *                * .                        * |
	 *      	     *  .                       *  |
	 *  	  	    *   .                      *   |
	 *  	  	   *    .                     *    |
	 * 	          *		.	                 *     |
	 *         5 *--------------------------* 7    |
	 *  	     |      .                   |      |
	 * 	     	 |      .                   |      |2
	 *           |	   0. . . . . . . . . . | .    *-----------> Y(j)
	 *           |     .                    |     *
	 *       	 |   .                      |   *
	 *  	     | .                        | *
	 *         4 *--------------------------*6
	 * 		  .
	        .
	      .
	     X(i)
	*/




	// Left Face
	mesh->faceToVertices[0]  = mesh->vertA;
	mesh->faceToVertices[1]  = mesh->vertE;
	mesh->faceToVertices[2]  = mesh->vertF;
	mesh->faceToVertices[3]  = mesh->vertB;

	// Front Face
	mesh->faceToVertices[4]  = mesh->vertE;
	mesh->faceToVertices[5]  = mesh->vertG;
	mesh->faceToVertices[6]  = mesh->vertH;
	mesh->faceToVertices[7]  = mesh->vertF;

	// Right Face
	mesh->faceToVertices[8]  = mesh->vertG;
	mesh->faceToVertices[9]  = mesh->vertC;
	mesh->faceToVertices[10] = mesh->vertD;
	mesh->faceToVertices[11] = mesh->vertH;

	// Back Face
	mesh->faceToVertices[12] = mesh->vertC;
	mesh->faceToVertices[13] = mesh->vertA;
	mesh->faceToVertices[14] = mesh->vertB;
	mesh->faceToVertices[15] = mesh->vertD;

	// Top Face
	mesh->faceToVertices[16] = mesh->vertB;
	mesh->faceToVertices[17] = mesh->vertF;
	mesh->faceToVertices[18] = mesh->vertH;
	mesh->faceToVertices[19] = mesh->vertD;

	// Bottom Face
	mesh->faceToVertices[20] = mesh->vertA;
	mesh->faceToVertices[21] = mesh->vertC;
	mesh->faceToVertices[22] = mesh->vertG;
	mesh->faceToVertices[23] = mesh->vertE;

//	meshFaceToVerticesConnectivityUnitTest(mesh);


}


/*************************************************************/
/* UNIT TEST                                                 */
/*************************************************************/

void meshCellToVerticesConnectivityUnitTest(mesh_t *mesh){

	/*                  ^ Z(k)
	 *                  .
	 *                  .
     *                  .
     *                 B*--------------------------*D
	 *				   *.                         *|
     *                * .                        * |
	 *      	     *  .                       *  |
	 *  	  	    *   .                      *   |
	 *  	  	   *    .                     *    |
	 * 	          *		.	                 *     |
	 *          F*--------------------------*H     |
	 *  	     |      .                   |      |
	 * 	     	 |      .                   |      |
	 *           |	   A. . . . . . . . . . | . . .*C-----------> Y(j)
	 *           |     .                    |     *
	 *       	 |   .                      |   *
	 *  	     | .                        | *
	 *          E*--------------------------*G
	 * 		  .
	        .
	      .
	     X(i)
	*/
	printf ( "\n!--------------------------------------------- !" );
	printf ( "\n! Cell To Vertex Connectivity                  !" );
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
	printf ( "\n! Cell To Vertex Connectivity                  !" );
	printf ( "\n!--------------------------------------------- !" );

	printf ( "\n%11s%5s%5s%5s%5s%5s%5s%5s%5s","! Vert ID :","A","B","C","D","E","F","G","H");
	printf ( "\n!--------------------------------------------- !\n" );
	for (dint cell = 0; cell < mesh->nicells; cell++) {

		printf ( "! Cell%3d :", cell );

		for (dint vert = 0; vert < mesh->ncellVertices; vert++) {

			dint vid = vert 	+ cell*mesh->ncellVertices;

			printf ( "%5d", mesh->cellToVertex[vid] );
		}

		printf ( "\n" );
	}

	printf ( "!--------------------------------------------- !\n" );


}

void meshCellToCellsConnectivityUnitTest(mesh_t *mesh){

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



	printf ( "\n!--------------------------------------------- !" );
	printf ( "\n! Cell To Cell Connectivity                    !" );
	printf ( "\n!--------------------------------------------- !" );

	printf ( "\n%11s%6s%7s%7s%6s%5s%8s","! Cell ID :","Left","Front","Right","Back","Top","Bottom");
	printf ( "\n!---------------------------------------------!" );


	for (dint cell = 0; cell < mesh->nicells; cell++) {

		printf ( "\n! Cell%3d :", cell );

		for (dint face = 0; face < mesh->ncellFaces; face++) {

			dint fid = face 	+ cell*mesh->ncellFaces;

			printf ( "%6d", mesh->cellToCells[fid] );
		}
	}

	printf ( "\n!---------------------------------------------!" );
}

void meshCellToFacesConnectivityUnitTest(mesh_t *mesh){

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


	printf ( "\n!--------------------------------------------- !" );
	printf ( "\n! Cell To Face Connectivity                    !" );
	printf ( "\n!--------------------------------------------- !" );

	printf ( "\n%11s%8s%8s%8s%8s%8s%8s","! Cell ID :","Left","Front","Right","Back","Top","Bottom");
	printf ( "\n!--------------------------------------------- !" );


	for (dint cell = 0; cell < mesh->nicells; cell++) {

		printf ( "\n! Cell%3d :", cell );

		for (dint face = 0; face < mesh->ncellFaces; face++) {

			dint fid = face 	+ cell*mesh->ncellFaces;

			if(mesh->cellToFaces[fid] == mesh->faceAEFB){ printf ( "%8s", "Left" );}
			if(mesh->cellToFaces[fid] == mesh->faceEGHF){ printf ( "%8s", "Front" );}
			if(mesh->cellToFaces[fid] == mesh->faceGCDH){ printf ( "%8s", "Right" );}
			if(mesh->cellToFaces[fid] == mesh->faceCABD){ printf ( "%8s", "Back" );}
			if(mesh->cellToFaces[fid] == mesh->faceFHDB){ printf ( "%8s", "Top" );}
			if(mesh->cellToFaces[fid] == mesh->faceACGE){ printf ( "%8s", "Bottom" );}
			if(mesh->cellToFaces[fid] == -1){ printf ( "%8s", "-1" );}

		}

	}
	printf ( "\n!--------------------------------------------- !\n" );

}


void meshFaceToVerticesConnectivityUnitTest(mesh_t *mesh){

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

	printf ( "\n*"
	"\n*                ^ Z(k)"
	"\n*                ."
	"\n*                ."
    "\n*                ."
    "\n*               1*--------------------------*3"
	"\n*               *                          *|"
    "\n*              * .                        * |"
	"\n*             *  .                       *  |"
	"\n*            *   .                      *   |"
	"\n*           *    .                     *    |"
	"\n*          *     .                    *     |"
	"\n*        5*--------------------------*7     |"
	"\n*         |      .                   |      |"
	"\n*         |      .                   |      |"
	"\n*         |     0. . . . . . . . . . | . . .*2-----------> Y(j)"
	"\n*         |     .                    |     *"
	"\n*         |   .                      |   *"
	"\n*         | .                        | *"
	"\n*        4*--------------------------*6"
	"\n*       ."
	"\n*     ."
	"\n*   ."
	"\n*  X(i)");

	printf ( "\n!--------------------------------------------- !" );
	printf ( "\n! Face To Vertices Connectivity                !" );
	printf ( "\n!--------------------------------------------- !" );
	printf ( "\n! Face ID :   Vertex1 Vertex2 Vertex3 Vertex4" );

	for (dint face = 0; face < mesh->ncellFaces; face++) {

		printf ( "\n! Face%3d : ", face );

		for (dint vertex = 0; vertex < mesh->nfaceVertices; vertex++) {

			dint vid = vertex 	+ face*mesh->nfaceVertices;

			printf ( " %6d", mesh->faceToVertices[vid] );
		}

	}
	printf ( "\n!--------------------------------------------- !\n" );


}

