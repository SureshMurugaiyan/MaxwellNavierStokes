#include "mesh.h"

void meshReadStatistics(mesh_t *mesh);

void meshReadStatistics(mesh_t *mesh){

	dint *meshStats;
	meshStats = (dint*) calloc(3, sizeof(dint));


	FILE *myFile;
	myFile = fopen("setUpCase/InputMesh/finiteVolumeMeshStatistics.txt", "r");


    if (myFile == NULL){
        printf("Error Reading File 'setUpCase/InputMesh/finiteVolumeMeshStatistics.txt' \n");
        printf("File Not Found \n");

        exit (0);
    }



    for (dint i = 0; i < 3; i++){
    	if (fscanf(myFile, "%d", &meshStats[i] ) > 0){

    	}else{
    		printf("cannot read data:Index exceeded:finiteVolumeMeshStatistics");
    	}
    }


    fclose(myFile);



	mesh->ncellsx 	 = meshStats[0];
	mesh->ncellsy 	 = meshStats[1];
	mesh->ncellsz 	 = meshStats[2];

	mesh->ncellVertices	 = 8;
	mesh->ncellFaces  	 = 6;
	mesh->nfaceVertices  = 4;
	mesh->ndimension 	 = 3;



	mesh->nicells  	 = mesh->ncellsx*mesh->ncellsy*mesh->ncellsz;


	// Calculating the number of boundary cells

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

	dfloat ncellsfaceAEFB = mesh->ncellsx*mesh->ncellsz;
	dfloat ncellsfaceGCDH = mesh->ncellsx*mesh->ncellsz;
	dfloat ncellsfaceCABD = mesh->ncellsy*mesh->ncellsz;
	dfloat ncellsfaceEGHF = mesh->ncellsy*mesh->ncellsz;
	dfloat ncellsfaceACGE = mesh->ncellsx*mesh->ncellsy;
	dfloat ncellsfaceFHDB = mesh->ncellsx*mesh->ncellsy;



	mesh->nbcells    = ncellsfaceAEFB + ncellsfaceGCDH + ncellsfaceCABD + ncellsfaceEGHF + ncellsfaceACGE + ncellsfaceFHDB;
	mesh->ncells  	 = mesh->nicells + mesh->nbcells;
	mesh->nivertices = mesh->nicells*mesh->ncellVertices;
	mesh->nvertices  = mesh->ncells*mesh->ncellVertices;

	mesh->dim_x = 0;
	mesh->dim_y = 1;
	mesh->dim_z = 2;

	free(meshStats);


//	printf ( "!--------------------------------- !\n" );
//	printf ( "! Mesh Statistics                  !\n" );
//	printf ( "!--------------------------------- !\n!" );
//	printf ( "%16s %5d \n!","ncellsx = ",mesh->ncellsx);
//	printf ( "%16s %5d \n!","ncellsy = ",mesh->ncellsy);
//	printf ( "%16s %5d \n!","ncellsz = ",mesh->ncellsz);
//
//	printf ( "%16s %5d \n!","ncellVertices = ",mesh->ncellVertices);
//	printf ( "%16s %5d \n!","ncellFaces = ",mesh->ncellFaces);
//	printf ( "%16s %5d \n!","nfaceVertices = ",mesh->nfaceVertices);
//	printf ( "%16s %5d \n!","ndimension = ",mesh->ndimension);
//
//
//	printf ( "%16s %5d \n!","nicells = ",mesh->nicells);
//	printf ( "%16s %5d \n!","nbcells = ",mesh->nbcells);
//	printf ( "%16s %5d \n!","ncells = ",mesh->ncells);
//	printf ( "%16s %5d \n!","nivertices = ",mesh->nivertices);
//	printf ( "%16s %5d \n!","nvertices = ",mesh->nvertices);
//	printf ( "%16s %5d \n!","dim_x = ",mesh->dim_x);
//	printf ( "%16s %5d \n!","dim_y = ",mesh->dim_y);
//	printf ( "%16s %5d \n","dim_z = ",mesh->dim_z);
//	printf ( "!--------------------------------- !\n" );

}
