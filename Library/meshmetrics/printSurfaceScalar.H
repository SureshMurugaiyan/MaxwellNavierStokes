#include "mesh.h"

void printSurfaceScalar(dfloat* Phi, mesh_t *mesh);

void printSurfaceScalar(dfloat* Phi, mesh_t *mesh){

	printMeshStencil(mesh);

	printf ( "!--------------------------------------------- !\n" );
	printf ( "!%10s                                    !\n",mesh->variableName);
	printf ( "!--------------------------------------------- !\n" );



	for (dint face = 0; face < mesh->ncellFaces; face++) {

		if (face == 0){printf ("!%12s : %1d\n!","Left Face(AEFB)",face);}
		if (face == 1){printf ("!%12s : %1d\n!","Front Face(EGHF)",face);}
		if (face == 2){printf ("!%12s : %1d\n!","Right Face(GCDH)",face);}
		if (face == 3){printf ("!%12s : %1d\n!","Back Face(CABD)",face);}
		if (face == 4){printf ("!%12s : %1d\n!","Top Face(FHDB)",face);}
		if (face == 5){printf ("!%12s : %1d\n!","Bottom Face(ACGE)",face);}

		printf ( "--------------------------------------------- !\n" );


		for( dint k = 0; k < mesh->ncellsz; k++){

			printf("! (Plane, K = %d)\n",k);

			for( dint i = 0; i < mesh->ncellsx; i++){

				printf("!");

				for( dint j = 0; j < mesh->ncellsy; j++){

					dint cell = computeLinearIndex(i,j,k, mesh->ncellsx, mesh->ncellsy,mesh->ncellsz);
					dint fid  = face+cell*mesh->ncellFaces;

					printf ( "%9.2e\t", Phi[fid]);


				}

				printf("\n");

			}
			printf ( "!--------------------------------------------- !\n" );
		}



	}

}
