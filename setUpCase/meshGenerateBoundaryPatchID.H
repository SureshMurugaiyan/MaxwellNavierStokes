#include "mesh.h"
#include "BoundaryNames.H"

void meshGenerateBoundaryPatchID(mesh_t *mesh);

void meshGenerateBoundaryPatchID(mesh_t *mesh){



	for (dint face = 0; face < mesh->nBoundaryFaces; face++) {

		dint fid  = mesh->boundaryFaces[face];

		dint Localface = computeLocalFace(fid,mesh);

		dfloat x = mesh->Cf[mesh->dim_x + fid*mesh->ndimension];
		dfloat y = mesh->Cf[mesh->dim_y + fid*mesh->ndimension];
		dfloat z = mesh->Cf[mesh->dim_z + fid*mesh->ndimension];

		if (Localface == mesh->faceAEFB) {mesh->boundaryFacePatchID[face] = leftBoundary;}
		if (Localface == mesh->faceEGHF) {mesh->boundaryFacePatchID[face] = frontBoundary;}
		if (Localface == mesh->faceGCDH) {mesh->boundaryFacePatchID[face] = rightBoundary;}
		if (Localface == mesh->faceCABD) {mesh->boundaryFacePatchID[face] = backBoundary;}
		if (Localface == mesh->faceFHDB) {mesh->boundaryFacePatchID[face] = topBoundary;}
		if (Localface == mesh->faceACGE) {mesh->boundaryFacePatchID[face] = bottomBoundary;}

	}


//	for (dint face = 0; face < mesh->nBoundaryFaces; face++) {
//
//		dint fid  = mesh->boundaryFaces[face];
//
//		dint Localface = computeLocalFace(fid,mesh);
//
//		dfloat x = mesh->Cf[mesh->dim_x + fid*mesh->ndimension];
//		dfloat y = mesh->Cf[mesh->dim_y + fid*mesh->ndimension];
//		dfloat z = mesh->Cf[mesh->dim_z + fid*mesh->ndimension];
//
//		if ((Localface == mesh->faceCABD) && (z>0.0)) {mesh->boundaryFacePatchID[face] = inlet;}
//
//		if ((Localface == mesh->faceEGHF) && (z>0.0)) {mesh->boundaryFacePatchID[face] = outlet;}
//
//	}


	// Template to assign boundaries based on Face centers
	//	for (dint face = 0; face < mesh->nBoundaryFaces; face++) {
	//
	//		dint fid = mesh->boundaryFaces[face];
	//
	//		dfloat x = mesh->Cf[mesh->dim_x + fid*mesh->ndimension];
	//		dfloat y = mesh->Cf[mesh->dim_y + fid*mesh->ndimension];
	//		dfloat z = mesh->Cf[mesh->dim_z + fid*mesh->ndimension];
	//
	//		if (magnitude(z-1.00) < mesh->SMALL) {mesh->boundaryFacePatchID[face] = topBoundary;   }
	//		if (magnitude(z+0.00) < mesh->SMALL) {mesh->boundaryFacePatchID[face] = bottomBoundary;}
	//		if (magnitude(y-0.00) < mesh->SMALL) {mesh->boundaryFacePatchID[face] = leftBoundary;  }
	//		if (magnitude(y-1.00) < mesh->SMALL) {mesh->boundaryFacePatchID[face] = rightBoundary; }
	//		if (magnitude(x-0.00) < mesh->SMALL) {mesh->boundaryFacePatchID[face] = backBoundary;  }
	//		if (magnitude(x-0.20) < mesh->SMALL) {mesh->boundaryFacePatchID[face] = frontBoundary; }
	//
	//	}








//	printMeshStencil(mesh);
//
//	printf ( "!--------------------------------------------- !\n" );
//	printf ( "! patchID                                      !\n" );
//	printf ( "!--------------------------------------------- !\n!" );
//	printf ( "%8s %20s\n","|FaceID|", "|patchID|");
//	printf ( "!--------------------------------------------- !\n!" );
//	for (dint face = 0; face < mesh->nBoundaryFaces; face++) {
//
//		printf ( "%4d %25s\n!", mesh->boundaryFaces[face],
//				getPatchID(mesh->boundaryFacePatchID[face]) );
//
//	}
//
//	printf ( "--------------------------------------------- !\n" );


}
