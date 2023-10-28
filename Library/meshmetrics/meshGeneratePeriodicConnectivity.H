#include "mesh.h"
#include "../setUpCase/BoundaryNames.H"

void meshGeneratePeriodicConnectivity(mesh_t *mesh);

void meshFaceToBoundaryCellConnectivity(mesh_t *mesh);
void meshPatchToPatchConnectivityPeriodic(mesh_t *mesh);
void meshFaceToNeighborCellConnectivityPeriodic(mesh_t *mesh);
void generateBoundaryFaceinterpolationFactorPeriodic(mesh_t *mesh);
void generateCellToCellDistanceVectorsPeriodic(mesh_t *mesh);
void meshCalculateNoOfInnerFacesIncludingPeriodicFaces(mesh_t *mesh);



void meshFaceToBoundaryCellConnectivityUnitTest(mesh_t *mesh);
void meshPatchToPatchConnectivityPeriodicUnitTest(mesh_t *mesh);
void meshFaceToNeighborCellConnectivityPeriodicUnitTest(mesh_t *mesh);
void generateBoundaryFaceinterpolationFactorPeriodicUnitTest(mesh_t *mesh);
void generateCellToCellDistanceVectorsPeriodicUnitTest(mesh_t *mesh);
void meshCalculateNoOfInnerFacesIncludingPeriodicFacesUnitTest(mesh_t *mesh);


void meshGeneratePeriodicConnectivity(mesh_t *mesh){

	meshFaceToBoundaryCellConnectivity(mesh);

	meshPatchToPatchConnectivityPeriodic(mesh);

	meshFaceToNeighborCellConnectivityPeriodic(mesh);

	//Modify the Metrics for Periodic connectivity
	generateBoundaryFaceinterpolationFactorPeriodic(mesh);

	generateCellToCellDistanceVectorsPeriodic(mesh);

	// Count of no of inner Faces including periodic faces
	meshCalculateNoOfInnerFacesIncludingPeriodicFaces(mesh);



}



void meshFaceToBoundaryCellConnectivity(mesh_t *mesh){


	for (dint cell = 0; cell < mesh->ncells; cell++) {
		for (dint face = 0; face < mesh->ncellFaces; face++) {

			dint fid 	 = face+cell*mesh->ncellFaces;
			dint cell_n  =  mesh->cellToCells[fid];

			mesh->faceToBoundaryCell[fid] = mesh->faceToNeighborCell[fid];

		}
	}

//	meshFaceToBoundaryCellConnectivityUnitTest(mesh);

}




void meshPatchToPatchConnectivityPeriodic(mesh_t *mesh){


	for (int patch = 0; patch < mesh->npatchs; patch++) {

		mesh->periodicNeighborPatch[patch] = -1; // Initialize with no periodic patch

		int patchID  = mesh->patchID[patch];
		int npatchID = mesh->periodicNeighborPatchID[patch];

		if(mesh->patchType[patch] == periodicPatch){

			for (int nbpatch = 0; nbpatch < mesh->npatchs; nbpatch++) {
				if (npatchID == mesh->patchID[nbpatch]){
					mesh->periodicNeighborPatch[patch] = nbpatch;
				}
			}
		}

	}


//	meshPatchToPatchConnectivityPeriodicUnitTest(mesh);

}




void meshFaceToNeighborCellConnectivityPeriodic(mesh_t *mesh){

	for (int patch = 0; patch < mesh->npatchs; patch++) {

		int nPatch = mesh->periodicNeighborPatch[patch];

		if ((mesh->patchType[patch] == periodicPatch)&& (mesh->npatchFaces[patch] !=mesh->npatchFaces[nPatch])){

			printf ( "\n" );
			printf ( "Fatal error!\n" );
			printf ( "  periodic boundaries should have equal number of cells!\n" );
			exit ( 1 );

		}
	}


	for (int patch = 0; patch < mesh->npatchs; patch++) {

		int npatch  = mesh->periodicNeighborPatch[patch];

		if(mesh->patchType[patch] == periodicPatch){

			int face    = mesh->patchStartFaceIndex[patch];
			int nface   = mesh->patchStartFaceIndex[npatch];

			for (int patchFace = 0; patchFace < mesh->npatchFaces[patch]; patchFace++){

				int fid    = mesh->patchToBoundaryFace[face];
				int nfid   = mesh->patchToBoundaryFace[nface];

				mesh->faceToNeighborCell[fid] = mesh->faceToOwnerCell[nfid];

				face++;
				nface++;
			}
		}
	}



//	meshFaceToNeighborCellConnectivityPeriodicUnitTest(mesh);

}

void generateBoundaryFaceinterpolationFactorPeriodic(mesh_t *mesh){

	for (int patch = 0; patch < mesh->npatchs; patch++) {

		int npatch  = mesh->periodicNeighborPatch[patch];

		if(mesh->patchType[patch] == periodicPatch){

			int face    = mesh->patchStartFaceIndex[patch];
			int nface   = mesh->patchStartFaceIndex[npatch];

			for (int patchFace = 0; patchFace < mesh->npatchFaces[patch]; patchFace++){

				int fid    = mesh->patchToBoundaryFace[face];
				int nfid   = mesh->patchToBoundaryFace[nface];

				dfloat PN1 = mesh->magdf[fid];
				dfloat PN2 = mesh->magdf[nfid];
				dfloat PN  = PN1+PN2;
				dfloat fN  = PN2;
				mesh->fx[fid] = fN/PN;


				face++;
				nface++;
			}
		}


	}

//	strcpy(mesh->variableName, "fx");
//	printSurfaceScalar(mesh->fx, mesh);


}



void generateCellToCellDistanceVectorsPeriodic(mesh_t *mesh){



	/*************************************************************/
	/* magnitude of d vector                                     */
	/*************************************************************/
	// copy to temp array

	dfloat *magdf_temp;
	createField(&magdf_temp,surfaceScalar, mesh);

	copyDataToFrom(magdf_temp,mesh->magdf,surfaceScalar,mesh);


	for (int patch = 0; patch < mesh->npatchs; patch++) {

		int npatch  = mesh->periodicNeighborPatch[patch];

		if(mesh->patchType[patch] == periodicPatch){

			int face    = mesh->patchStartFaceIndex[patch];
			int nface   = mesh->patchStartFaceIndex[npatch];

			for (int patchFace = 0; patchFace < mesh->npatchFaces[patch]; patchFace++){

				int fid    = mesh->patchToBoundaryFace[face];
				int nfid   = mesh->patchToBoundaryFace[nface];


				mesh->magd[fid]  =  magdf_temp[fid] +  magdf_temp[nfid] ;

				mesh->rmagd[fid] = 1.0/mesh->magd[fid];

				mesh->delta[fid]   = mesh->magd[fid];
				mesh->rdelta[fid]  = mesh->rmagd[fid];

				face++;
				nface++;
			}
		}

	}

	free(magdf_temp);



	/*************************************************************/
	/*  d vector                                                 */
	/*************************************************************/

	for (int patch = 0; patch < mesh->npatchs; patch++) {

		int npatch  = mesh->periodicNeighborPatch[patch];

		if(mesh->patchType[patch] == periodicPatch){

			int face    = mesh->patchStartFaceIndex[patch];

			for (int patchFace = 0; patchFace < mesh->npatchFaces[patch]; patchFace++){

				int fid    = mesh->patchToBoundaryFace[face];

				dfloat *d_cell;
				createField(&d_cell,vector, mesh);

				for (int dim = 0; dim < mesh->ndimension; dim++) {
					int fid_dim   = dim + fid*mesh->ndimension;

					d_cell[dim] = mesh->df[fid_dim] ;
				}

				dfloat magD = magVector(d_cell);

				free(d_cell);

				for (int dim = 0; dim < mesh->ndimension; dim++) {
					int fid_dim   = dim + fid*mesh->ndimension;

					mesh->d[fid_dim] = (mesh->df[fid_dim]/magD)*mesh->magd[fid] ;

				}


				face++;
			}
		}

	}

//	strcpy(mesh->variableName, "d");
//	printSurfaceVector(mesh->d, mesh);
//
//	strcpy(mesh->variableName, "magd");
//	printSurfaceScalar(mesh->magd, mesh);
//
//	strcpy(mesh->variableName, "rmagd");
//	printSurfaceScalar(mesh->rmagd, mesh);
//
//	strcpy(mesh->variableName, "delta");
//	printSurfaceScalar(mesh->delta, mesh);
//
//	strcpy(mesh->variableName, "rdelta");
//	printSurfaceScalar(mesh->rdelta, mesh);

}




void meshCalculateNoOfInnerFacesIncludingPeriodicFaces(mesh_t *mesh){

	mesh->nInnerAndPeriodicFaces = 0;
	mesh->nBoundaryFaces = 0;

	for (int cell = 0; cell < mesh->nicells; cell++) {
		for (int face = 0; face < mesh->ncellFaces; face++) {

			int fid = face+cell*mesh->ncellFaces;
			int cell_n = mesh->faceToNeighborCell[fid];

			if (cell_n < mesh->nicells) {

				mesh->nInnerAndPeriodicFaces = mesh->nInnerAndPeriodicFaces+1;
			}else{

				mesh->nBoundaryFaces = mesh->nBoundaryFaces+1;
			}


		}
	}



//	printf ( "!----------------------------------------------- !\n" );
//	printf ( "! Number of Inner Faces Including Periodic Faces !\n" );
//	printf ( "!----------------------------------------------- !\n!" );
//
//	printf(" No of Inner Faces and Periodic Faces, "
//			"N = %d\t \n!",mesh->nInnerAndPeriodicFaces );
//	printf ( "----------------------------------------------- !\n" );
//
//	printf ( "!--------------------------------------------- !\n" );
//	printf ( "! Number of Boundary Faces                     !\n" );
//	printf ( "!--------------------------------------------- !\n!" );
//
//	printf(" No of Boundary Faces,     N = %d\t \n!", mesh->nBoundaryFaces);
//	printf ( "--------------------------------------------- !\n" );

}



/*************************************************************/
/* UNIT TEST                                                 */
/*************************************************************/

void meshPatchToPatchConnectivityPeriodicUnitTest(mesh_t *mesh){

	printf ( "!--------------------------------------------- !\n" );
	printf ( "! neighborPatchID of Each Patch                !\n" );
	printf ( "!--------------------------------------------- !\n!" );

	printf ( " %16s  %8s \n!","Patch ID","neighborPatch ID" );
	printf ( "--------------------------------------------- !\n!" );
	for (int patch = 0; patch < mesh->npatchs; patch++) {
		int nPatch = mesh->periodicNeighborPatch[patch];

		if(nPatch > -1){
			printf ( " %18s(%d)  %18s(%d) \n!", getPatchID(mesh->patchID[patch]) ,mesh->patchID[patch],
					getPatchID (mesh->patchID[nPatch]),mesh->patchID[nPatch]);
		}else{

			printf ( " %18s(%d)  %18s(%d) \n!", getPatchID(mesh->patchID[patch]) ,mesh->patchID[patch],
					"none",nPatch);

		}
	}

	printf ( "--------------------------------------------- !\n" );



}



void meshFaceToNeighborCellConnectivityPeriodicUnitTest(mesh_t *mesh){


	printf ( "!--------------------------------------------- !\n" );
	printf ( "! Face to NeighborCell   Connectivity          !\n" );
	printf ( "!--------------------------------------------- !\n" );

	printf ( "%11s%6s%6s%6s%6s\n","! Cell ID :","West","South","East","North");
	printf ( "!--------------------------------------------- !\n!" );


	for (int cell = 0; cell < mesh->nicells; cell++) {

		printf ( " Cell%3d :", cell );

		for (int face = 0; face < mesh->ncellFaces; face++) {

			int fid = face 	+ cell*mesh->ncellFaces;

			printf ( "%6d", mesh->faceToNeighborCell[fid] );
		}

		printf ( "\n!" );
	}

	printf ( "--------------------------------------------- !\n" );

}


void meshFaceToBoundaryCellConnectivityUnitTest(mesh_t *mesh){

	printMeshStencil(mesh);

	printf ( "!--------------------------------------------- !\n" );
	printf ( "! Face to BoundaryCell Connectivity            !\n" );
	printf ( "!--------------------------------------------- !\n" );

	printf ( "%11s%6s%7s%7s%6s%5s%8s\n","! Cell ID :","Left","Front","Right","Back","Top","Bottom");
	printf ( "!--------------------------------------------- !\n!" );


	for (dint cell = 0; cell < mesh->nicells; cell++) {

		printf ( " Cell%3d :", cell );

		for (dint face = 0; face < mesh->ncellFaces; face++) {

			dint fid = face 	+ cell*mesh->ncellFaces;

			printf ( "%6d", mesh->faceToBoundaryCell[fid] );
		}

		printf ( "\n!" );
	}

	printf ( "--------------------------------------------- !\n" );

}

