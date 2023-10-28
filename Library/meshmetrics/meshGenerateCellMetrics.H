#include "mesh.h"

void meshGenerateCellMetrics(mesh_t *mesh);

void generateCellCenters(mesh_t *mesh);
void generateFaceCenters(mesh_t *mesh);
void generateFaceAreas(mesh_t *mesh);
void generateFaceNormaVectors(mesh_t *mesh);
void generateCellVolumes(mesh_t *mesh);
void generateCellToCellDistanceVectors(mesh_t *mesh);
void generateCellToFaceDistanceVectors(mesh_t *mesh);
void generateDeltaCoefficients(mesh_t *mesh);
void generateCellFaceInterpolationFactor(mesh_t *mesh);
void generateCellCentersBoundaryCells(mesh_t *mesh);


void meshGenerateCellMetrics(mesh_t *mesh){

	generateCellCenters(mesh);

	generateFaceCenters(mesh);

	generateFaceAreas(mesh);

	generateFaceNormaVectors(mesh);

	generateCellVolumes(mesh);

	generateCellToCellDistanceVectors(mesh);

	generateCellToFaceDistanceVectors(mesh);

	generateDeltaCoefficients(mesh);

	generateCellFaceInterpolationFactor(mesh);

	generateCellCentersBoundaryCells(mesh);



}

void generateCellCenters(mesh_t *mesh){

	dfloat *sum;
	sum = (dfloat*) calloc(mesh->ndimension, sizeof(dfloat));

	for( dint cell = 0; cell < mesh->nicells;cell++){

		// Clear container

		for (dint i = 0; i < mesh->ndimension; i++) {sum[i] = 0;}


		for(dint vertex=0;vertex<mesh->ncellVertices;vertex++){
			for(dint dimension=0;dimension<mesh->ndimension;dimension++){

				dint vid_s 		= vertex + cell*mesh->ncellVertices;
				dint vid 		= mesh->cellToVertex[vid_s];
				dint vid_dim	= dimension + vid*mesh->ndimension;


				sum[dimension] = sum[dimension] + mesh->vertices[vid_dim];

			}
		}

		for(dint dimension=0;dimension<mesh->ndimension;dimension++){
			dint cell_dim		= dimension + cell*mesh->ndimension;
			mesh->C[cell_dim] 	= sum[dimension]/(dfloat)mesh->ncellVertices;
		}

	}

	free(sum);


//	strcpy(mesh->variableName, "C");
//	printVolumeVector(mesh->C, mesh);

}



void generateFaceCenters(mesh_t *mesh){

	dfloat *sum;
	sum = (dfloat*) calloc(mesh->ndimension, sizeof(dfloat));

	for( dint cell = 0; cell < mesh->nicells;cell++){
		for (dint face = 0; face < mesh->ncellFaces;face++){

			// Clear container
			for (dint i = 0; i < mesh->ndimension; i++) {sum[i] = 0;}

			for (dint facevertex = 0; facevertex < mesh->nfaceVertices; facevertex++){

				dint vid_ss = facevertex+face*mesh->nfaceVertices;
				dint vertex = mesh->faceToVertices[vid_ss];
				dint vid    = vertex + cell*mesh->ncellVertices;

				for(dint dimension=0;dimension<mesh->ndimension;dimension++){

					dint vid_dim   = dimension + vid*mesh->ndimension;
					sum[dimension] = sum[dimension] + mesh->vertices[vid_dim];

				}
			}

			dint fid  = face + cell*mesh->ncellFaces;

			for(dint dimension=0;dimension<mesh->ndimension;dimension++){

				dint fid_dim  = dimension + fid*mesh->ndimension;

				mesh->Cf[fid_dim] = sum[dimension]/(dfloat)mesh->nfaceVertices;

			}

		}
	}

	free(sum);


//	strcpy(mesh->variableName, "Cf");
//	printSurfaceVector(mesh->Cf, mesh);

}

void generateFaceAreas(mesh_t *mesh){

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

	dfloat *R02;
	dfloat *R13;

	R02 = (dfloat*) calloc(mesh->ndimension, sizeof(dfloat));
	R13 = (dfloat*) calloc(mesh->ndimension, sizeof(dfloat));


	for( dint cell = 0; cell < mesh->nicells;cell++) {
		for (dint face = 0; face < mesh->ncellFaces;face++) {
			for (dint dim = 0; dim < mesh->ndimension; dim++) {

				dint fid = face + cell*mesh->ncellFaces;

				dint vert0 = 0;
				dint vert1 = 1;
				dint vert2 = 2;
				dint vert3 = 3;

				dint vid_0_s = mesh->faceToVertices[vert0 + face*mesh->nfaceVertices];
				dint vid_1_s = mesh->faceToVertices[vert1 + face*mesh->nfaceVertices];
				dint vid_2_s = mesh->faceToVertices[vert2 + face*mesh->nfaceVertices];
				dint vid_3_s = mesh->faceToVertices[vert3 + face*mesh->nfaceVertices];

				dint vid_0 = vid_0_s + cell*mesh->ncellVertices;
				dint vid_1 = vid_1_s + cell*mesh->ncellVertices;
				dint vid_2 = vid_2_s + cell*mesh->ncellVertices;
				dint vid_3 = vid_3_s + cell*mesh->ncellVertices;

				dint vid_0_dim  = dim + vid_0*mesh->ndimension;
				dint vid_1_dim  = dim + vid_1*mesh->ndimension;
				dint vid_2_dim  = dim + vid_2*mesh->ndimension;
				dint vid_3_dim  = dim + vid_3*mesh->ndimension;

				R02[dim] = mesh->vertices[vid_2_dim] - mesh->vertices[vid_0_dim];
				R13[dim] = mesh->vertices[vid_3_dim] - mesh->vertices[vid_1_dim];

			}

			dint fid  = face + cell*mesh->ncellFaces;

			dint fid_x = mesh->dim_x + fid*mesh->ndimension;
			dint fid_y = mesh->dim_y + fid*mesh->ndimension;
			dint fid_z = mesh->dim_z + fid*mesh->ndimension;

			dfloat a1 = R02[0];
			dfloat a2 = R02[1];
			dfloat a3 = R02[2];

			dfloat b1 = R13[0];
			dfloat b2 = R13[1];
			dfloat b3 = R13[2];

			mesh->Sf[fid_x] = 0.5*(a2*b3 - a3*b2);
			mesh->Sf[fid_y] = 0.5*(a3*b1 - a1*b3);
			mesh->Sf[fid_z] = 0.5*(a1*b2 - a2*b1);

		}
	}



	free(R02);
	free(R13);


	dfloat *area;
	area = (dfloat*) calloc(mesh->ndimension, sizeof(dfloat));

	// Compute the magnitude of the face areas
	for( dint cell = 0; cell < mesh->nicells;cell++){
		for (dint face = 0; face < mesh->ncellFaces; face++) {

			dint fid = face + cell*mesh->ncellFaces;

			for (dint dim = 0; dim < mesh->ndimension; dim++) {

				dint fid_dim = dim + fid*mesh->ndimension;
				area[dim] = mesh->Sf[fid_dim];
			}

			mesh->magSf[fid] = sqrt(pow(area[0],2) + pow(area[1],2) + pow(area[2],2));

		}
	}

	free(area);



//	strcpy(mesh->variableName, "Sf");
//	printSurfaceVector(mesh->Sf, mesh);
//
//	strcpy(mesh->variableName, "magSf");
//	printSurfaceScalar(mesh->magSf, mesh);

}


void generateFaceNormaVectors(mesh_t *mesh){

	for( dint cell = 0; cell < mesh->nicells;cell++){
		for (dint face = 0; face < mesh->ncellFaces; ++face) {

			dint fid = face + cell*mesh->ncellFaces;

			dint fid_x = 0 + fid*mesh->ndimension;
			dint fid_y = 1 + fid*mesh->ndimension;
			dint fid_z = 2 + fid*mesh->ndimension;

			dfloat a = mesh->Sf[fid_x];
			dfloat b = mesh->Sf[fid_y];
			dfloat c = mesh->Sf[fid_z];

			dfloat magVector = sqrt(a*a +b*b +c*c);


			for (dint dimension = 0; dimension < mesh->ndimension; dimension++) {

				dint fid_dim = dimension + fid*mesh->ndimension;

				magVector = stabilise(magVector,mesh);

				mesh->Sn[fid_dim] = mesh->Sf[fid_dim]/magVector ;

			}

		}
	}

//	strcpy(mesh->variableName, "Sn");
//	printSurfaceVector(mesh->Sn, mesh);


}


void generateCellVolumes(mesh_t *mesh){


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


	dfloat *S_faceGCDH; // right face
	dfloat *S_faceCABD; // back face
	dfloat *S_faceFHDB; // top face
	dfloat *sum;
	dfloat *rD_minus_rE;

	S_faceGCDH  = (dfloat*) calloc(mesh->ndimension, sizeof(dfloat));
	S_faceCABD  = (dfloat*) calloc(mesh->ndimension, sizeof(dfloat));
	S_faceFHDB  = (dfloat*) calloc(mesh->ndimension, sizeof(dfloat));
	sum         = (dfloat*) calloc(mesh->ndimension, sizeof(dfloat));
	rD_minus_rE = (dfloat*) calloc(mesh->ndimension, sizeof(dfloat));



	for( dint cell = 0; cell < mesh->nicells;cell++) {


		dint fid_faceGCDH = mesh->faceGCDH + cell*mesh->ncellFaces; // right face
		dint fid_faceCABD = mesh->faceCABD + cell*mesh->ncellFaces; // back face
		dint fid_faceFHDB = mesh->faceFHDB + cell*mesh->ncellFaces; // top face

		for (dint dim = 0; dim < mesh->ndimension; dim++) {

			dint fid_faceGCDH_dim =  dim + fid_faceGCDH*mesh->ndimension;
			dint fid_faceCABD_dim =  dim + fid_faceCABD*mesh->ndimension;
			dint fid_faceFHDB_dim =  dim + fid_faceFHDB*mesh->ndimension;

			S_faceGCDH[dim] = mesh->Sf[fid_faceGCDH_dim];
			S_faceCABD[dim] = mesh->Sf[fid_faceCABD_dim];
			S_faceFHDB[dim] = mesh->Sf[fid_faceFHDB_dim];

		}

		for (dint dim = 0; dim < mesh->ndimension; dim++) {

			sum[dim] = S_faceGCDH[dim] + S_faceCABD[dim] + S_faceFHDB[dim];

		}


		for (dint dim = 0; dim < mesh->ndimension; dim++) {

			dint vid_vertD = mesh->vertD + cell*mesh->ncellVertices;
			dint vid_vertE = mesh->vertE + cell*mesh->ncellVertices;

			dint vertD_dim = dim + vid_vertD*mesh->ndimension;
			dint vertE_dim = dim + vid_vertE*mesh->ndimension;

			rD_minus_rE[dim] = mesh->vertices[vertD_dim] - mesh->vertices[vertE_dim];
		}

		mesh->V[cell]  =  (1.0/3.0)*(sum[0]*rD_minus_rE[0]+ sum[1]*rD_minus_rE[1] + sum[2]*rD_minus_rE[2]);

		mesh->rV[cell] = 1.0/mesh->V[cell] ; // reciprocal of Volumes
	}

	free(S_faceGCDH);
	free(S_faceCABD);
	free(S_faceFHDB);
	free(sum);
	free(rD_minus_rE);


//	strcpy(mesh->variableName, "Vol");
//	printVolumeScalar(mesh->V, mesh);
//
//	strcpy(mesh->variableName, "rV");
//	printVolumeScalar(mesh->rV, mesh);

}




void generateCellToCellDistanceVectors(mesh_t *mesh){


	for (dint cell = 0; cell < mesh->nicells; ++cell) {
		for (dint face = 0; face < mesh->ncellFaces; ++face) {

			dint fid	=  face  + cell*mesh->ncellFaces;
			dint cell_n	=  mesh->faceToNeighborCell[fid];

			if(mesh->cellToFaces[fid] >(-1)){

				for (dint dim = 0; dim < mesh->ndimension; dim++) {

					dint cell_dim  = dim + cell*mesh->ndimension;
					dint cell_ndim = dim + cell_n*mesh->ndimension;
					dint fid_dim   = dim + fid*mesh->ndimension;

					mesh->d[fid_dim] = mesh->C[cell_ndim]-mesh->C[cell_dim];

				}
			}

		}
	}




	// Compute the magnitude of the distance vector

	dfloat *dl;
	dl = (dfloat*) calloc(mesh->ndimension, sizeof(dfloat));

	for (dint cell = 0; cell < mesh->nicells; ++cell) {
		for (dint face = 0; face < mesh->ncellFaces; ++face) {

			dint fid =  face  + cell*mesh->ncellFaces;

			for (dint dim = 0; dim < mesh->ndimension; dim++) {
				dint fid_dim   = dim + fid*mesh->ndimension;
				dl[dim] = mesh->d[fid_dim];
			}

			if(mesh->cellToFaces[fid] >(-1)){

				mesh->magd[fid] = sqrt(pow(dl[0],2) + pow(dl[1],2) + pow(dl[2],2));

				mesh->rmagd[fid] = 1.0/mesh->magd[fid];
			}

		}
	}

	free(dl);



//	strcpy(mesh->variableName, "d");
//	printSurfaceVector(mesh->d, mesh);
//
//	strcpy(mesh->variableName, "magd");
//	printSurfaceScalar(mesh->magd, mesh);
//
//	strcpy(mesh->variableName, "rmagd");
//	printSurfaceScalar(mesh->rmagd, mesh);


}


void generateCellToFaceDistanceVectors(mesh_t *mesh){


	for (dint cell = 0; cell < mesh->nicells; ++cell) {
		for (dint face = 0; face < mesh->ncellFaces; ++face) {

			dint fid	=  face  + cell*mesh->ncellFaces;
			dint cell_n	=  mesh->faceToNeighborCell[fid];

			for (dint dim = 0; dim < mesh->ndimension; dim++) {

				dint cell_dim  = dim + cell*mesh->ndimension;
				dint cell_ndim = dim + cell_n*mesh->ndimension;
				dint fid_dim   = dim + fid*mesh->ndimension;

				mesh->df[fid_dim] = mesh->Cf[fid_dim]-mesh->C[cell_dim];

			}

		}
	}



	// Compute the magnitude of the distance vector

	dfloat *dl;
	dl = (dfloat*) calloc(mesh->ndimension, sizeof(dfloat));

	for (dint cell = 0; cell < mesh->nicells; ++cell) {
		for (dint face = 0; face < mesh->ncellFaces; ++face) {

			dint fid		=  face  + cell*mesh->ncellFaces;

			for (dint dim = 0; dim < mesh->ndimension; dim++) {
				dint fid_dim   = dim + fid*mesh->ndimension;
				dl[dim] = mesh->df[fid_dim];
			}

			mesh->magdf[fid] = sqrt(pow(dl[0],2) + pow(dl[1],2) + pow(dl[2],2));

			mesh->rmagdf[fid] = 1.0/mesh->magdf[fid];
		}
	}

	free(dl);


//	strcpy(mesh->variableName, "df");
//	printSurfaceVector(mesh->df, mesh);
//
//	strcpy(mesh->variableName, "magdf");
//	printSurfaceScalar(mesh->magdf, mesh);
//
//	strcpy(mesh->variableName, "rmagdf");
//	printSurfaceScalar(mesh->rmagdf, mesh);

}

void generateDeltaCoefficients(mesh_t *mesh){


	for (dint cell = 0; cell < mesh->nicells; ++cell) {
		for (dint face = 0; face < mesh->ncellFaces; ++face) {

			dint fid =  face  + cell*mesh->ncellFaces;


			if(mesh->cellToFaces[fid] >(-1)){

				mesh->delta[fid]  = mesh->magd[fid];
				mesh->rdelta[fid] = mesh->rmagd[fid];


			} else{

				mesh->delta[fid]  = mesh->magdf[fid];
				mesh->rdelta[fid] = mesh->rmagdf[fid];

			}

		}
	}
//
//	strcpy(mesh->variableName, "delta");
//	printSurfaceScalar(mesh->delta, mesh);
//
//	strcpy(mesh->variableName, "rdelta");
//	printSurfaceScalar(mesh->rdelta, mesh);



}


void generateCellFaceInterpolationFactor(mesh_t *mesh){


	for (dint cell = 0; cell < mesh->nicells; ++cell) {
		for (dint face = 0; face < mesh->ncellFaces; ++face) {

			dint fid	=  face  + cell*mesh->ncellFaces;
			dint cell_n	=  mesh->faceToNeighborCell[fid];

			if(mesh->cellToFaces[fid] >(-1)){
				mesh->fx[fid] = mesh->magdf[fid]/mesh->magd[fid];
			}


		}
	}

//	strcpy(mesh->variableName, "fx");
//	printSurfaceScalar(mesh->fx, mesh);


}

void generateCellCentersBoundaryCells(mesh_t *mesh){


	for (dint cell = 0; cell < mesh->nicells; ++cell) {
		for (dint face = 0; face < mesh->ncellFaces; ++face) {

			dint fid	=  face  + cell*mesh->ncellFaces;
			dint cell_n	=  mesh->faceToNeighborCell[fid];

			if(mesh->cellToFaces[fid] ==(-1)){

				for(dint dimension=0;dimension<mesh->ndimension;dimension++){

					dint cell_n_dim  = dimension + cell_n*mesh->ndimension;
					dint fid_dim     = dimension + fid*mesh->ndimension;

					mesh->C[cell_n_dim] = mesh->Cf[fid_dim];
				}

			}


		}
	}


//	strcpy(mesh->variableName, "C");
//	printVolumeVector(mesh->C, mesh);

}
