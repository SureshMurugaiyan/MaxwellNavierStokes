#include "mesh.h"

void printVolumeVector(dfloat* Phi, mesh_t *mesh);

void printVolumeVector(dfloat* Phi, mesh_t *mesh){

	printMeshStencil(mesh);

	printf ( "\n!--------------------------------------------- !" );
	printf ( "\n!%3s                                           !",mesh->variableName);
	printf ( "\n!--------------------------------------------- !" );



	printf ( "\n!--------------------------------------------- !" );
	printf ( "\n! Inner cells                                  !" );
	printf ( "\n!--------------------------------------------- !" );



	for( dint k = (mesh->ncellsz-1); k >=0; k--){

		printf("\n! (Plane, K = %d)\n",k);

		for( dint i = 0; i < mesh->ncellsx; i++){

			printf("!");

			for( dint j = 0; j < mesh->ncellsy; j++){

				dint cell = computeLinearIndex(i,j,k, mesh->ncellsx, mesh->ncellsy,mesh->ncellsz);

				dint id_x = mesh->dim_x + cell*mesh->ndimension;
				dint id_y = mesh->dim_y + cell*mesh->ndimension;
				dint id_z = mesh->dim_z + cell*mesh->ndimension;

				printf ( " (%9.2e,%9.2e,%9.2e)", Phi[id_x], Phi[id_y], Phi[id_z] );


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

			dint id_x = mesh->dim_x + bcell*mesh->ndimension;
			dint id_y = mesh->dim_y + bcell*mesh->ndimension;
			dint id_z = mesh->dim_z + bcell*mesh->ndimension;

			printf ( " (%9.2e,%9.2e,%9.2e)", Phi[id_x], Phi[id_y], Phi[id_z] );


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
			dint bcell = mesh->faceToNeighborCell[fid];

			dint id_x = mesh->dim_x + bcell*mesh->ndimension;
			dint id_y = mesh->dim_y + bcell*mesh->ndimension;
			dint id_z = mesh->dim_z + bcell*mesh->ndimension;

			printf ( " (%9.2e,%9.2e,%9.2e)", Phi[id_x], Phi[id_y], Phi[id_z] );


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

			dint id_x = mesh->dim_x + bcell*mesh->ndimension;
			dint id_y = mesh->dim_y + bcell*mesh->ndimension;
			dint id_z = mesh->dim_z + bcell*mesh->ndimension;

			printf ( " (%9.2e,%9.2e,%9.2e)", Phi[id_x], Phi[id_y], Phi[id_z] );


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

			dint id_x = mesh->dim_x + bcell*mesh->ndimension;
			dint id_y = mesh->dim_y + bcell*mesh->ndimension;
			dint id_z = mesh->dim_z + bcell*mesh->ndimension;

			printf ( " (%9.2e,%9.2e,%9.2e)", Phi[id_x], Phi[id_y], Phi[id_z] );


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

			dint id_x = mesh->dim_x + bcell*mesh->ndimension;
			dint id_y = mesh->dim_y + bcell*mesh->ndimension;
			dint id_z = mesh->dim_z + bcell*mesh->ndimension;

			printf ( " (%9.2e,%9.2e,%9.2e)", Phi[id_x], Phi[id_y], Phi[id_z] );


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

			dint id_x = mesh->dim_x + bcell*mesh->ndimension;
			dint id_y = mesh->dim_y + bcell*mesh->ndimension;
			dint id_z = mesh->dim_z + bcell*mesh->ndimension;

			printf ( " (%9.2e,%9.2e,%9.2e)", Phi[id_x], Phi[id_y], Phi[id_z] );


		}
		printf("\n");
	}

	printf ( "!--------------------------------------------- !\n" );



}
