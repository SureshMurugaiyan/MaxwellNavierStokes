
#include "mesh.h"

void printVolumeScalar(dfloat* Phi, mesh_t *mesh);

void printVolumeScalar(dfloat* Phi, mesh_t *mesh){

	printMeshStencil(mesh);

	printf ( "!--------------------------------------------- !\n" );
	printf ( "!%3s                                           !\n",mesh->variableName);
	printf ( "!--------------------------------------------- !\n" );


	printf ( "!--------------------------------------------- !\n" );
	printf ( "! Inner cells                                  !\n" );
	printf ( "!--------------------------------------------- !\n" );



	for( dint k = (mesh->ncellsz-1); k >=0; k--){

		printf("! (Plane, K = %d)\n",k);

		for( dint i = 0; i < mesh->ncellsx; i++){

			printf("!");

			for( dint j = 0; j < mesh->ncellsy; j++){

				dint cell = computeLinearIndex(i,j,k, mesh->ncellsx, mesh->ncellsy,mesh->ncellsz);

				printf ( "%9.3e\t", Phi[cell]);

			}

			printf("\n");

		}
		printf ( "!--------------------------------------------- !\n" );
	}



	printf ( "!--------------------------------------------- !\n" );
	printf ( "! Bottom Boundary cells                        !\n" );
	printf ( "!--------------------------------------------- !\n" );


	for( dint i = 0; i < mesh->ncellsx; i++){

		printf("!");

		for( dint j = 0; j < mesh->ncellsy; j++){

			dint k = 0;

			dint cell  = computeLinearIndex(i,j,k, mesh->ncellsx, mesh->ncellsy, mesh->ncellsz);
			dint fid   = mesh->faceACGE + cell*mesh->ncellFaces;
			dint bcell = mesh->faceToNeighborCell[fid] ;

			printf ( "%9.3e\t", Phi[bcell]);


		}
		printf("\n");
	}

	printf ( "!--------------------------------------------- !\n" );

	printf ( "!--------------------------------------------- !\n" );
	printf ( "! Top Boundary cells                           !\n" );
	printf ( "!--------------------------------------------- !\n" );


	for( dint i = 0; i < mesh->ncellsx; i++){

		printf("!");

		for( dint j = 0; j < mesh->ncellsy; j++){

			dint k = mesh->ncellsz-1;

			dint cell  = computeLinearIndex(i,j,k, mesh->ncellsx, mesh->ncellsy, mesh->ncellsz);
			dint fid   = mesh->faceFHDB + cell*mesh->ncellFaces;
			dint bcell = mesh->faceToNeighborCell[fid] ;

			printf ( "%9.3e\t", Phi[bcell]);


		}
		printf("\n");
	}

	printf ( "!--------------------------------------------- !\n" );

	printf ( "!--------------------------------------------- !\n" );
	printf ( "! Left Boundary cells                          !\n" );
	printf ( "!--------------------------------------------- !\n" );


	for( dint k = (mesh->ncellsz-1); k >=0; k--){

		printf("!");

		for( dint i = 0; i < mesh->ncellsx; i++){

			dint j = 0;

			dint cell  = computeLinearIndex(i,j,k, mesh->ncellsx, mesh->ncellsy, mesh->ncellsz);
			dint fid   = mesh->faceAEFB + cell*mesh->ncellFaces;
			dint bcell = mesh->faceToNeighborCell[fid] ;

			printf ( "%9.3e\t", Phi[bcell]);


		}
		printf("\n");
	}

	printf ( "!--------------------------------------------- !\n" );

	printf ( "!--------------------------------------------- !\n" );
	printf ( "! Right Boundary cells                          !\n" );
	printf ( "!--------------------------------------------- !\n" );


	for( dint k = (mesh->ncellsz-1); k >=0; k--){

		printf("!");

		for( dint i = 0; i < mesh->ncellsx; i++){

			dint j = mesh->ncellsy-1;

			dint cell  = computeLinearIndex(i,j,k, mesh->ncellsx, mesh->ncellsy, mesh->ncellsz);
			dint fid   = mesh->faceGCDH + cell*mesh->ncellFaces;
			dint bcell = mesh->faceToNeighborCell[fid] ;

			printf ( "%9.3e\t", Phi[bcell]);


		}
		printf("\n");
	}

	printf ( "!--------------------------------------------- !\n" );

	printf ( "!--------------------------------------------- !\n" );
	printf ( "! Back Boundary cells                          !\n" );
	printf ( "!--------------------------------------------- !\n" );


	for( dint k = (mesh->ncellsz-1); k >=0; k--){

		printf("!");

		for( dint j = 0; j < mesh->ncellsy; j++){

			dint i = 0;

			dint cell  = computeLinearIndex(i,j,k, mesh->ncellsx, mesh->ncellsy, mesh->ncellsz);
			dint fid   = mesh->faceCABD + cell*mesh->ncellFaces;
			dint bcell = mesh->faceToNeighborCell[fid] ;

			printf ( "%9.3e\t", Phi[bcell]);


		}
		printf("\n");
	}

	printf ( "!--------------------------------------------- !\n" );

	printf ( "!--------------------------------------------- !\n" );
	printf ( "! Front Boundary cells                         !\n" );
	printf ( "!--------------------------------------------- !\n" );


	for( dint k = (mesh->ncellsz-1); k >=0; k--){

		printf("!");

		for( dint j = 0; j < mesh->ncellsy; j++){

			dint i = mesh->ncellsx-1;

			dint cell  = computeLinearIndex(i,j,k, mesh->ncellsx, mesh->ncellsy, mesh->ncellsz);
			dint fid   = mesh->faceEGHF + cell*mesh->ncellFaces;
			dint bcell = mesh->faceToNeighborCell[fid] ;

			printf ( "%9.3e\t", Phi[bcell]);


		}
		printf("\n");
	}

	printf ( "!--------------------------------------------- !\n" );



}
