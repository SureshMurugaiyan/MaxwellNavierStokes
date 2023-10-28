#include "mesh.h"

//#include "../Library/fvs/extrapolateToBoundaryCellVector.H"

void extrapolateToBoundaryCellVector(dfloat* phiCell,mesh_t *mesh);

void extrapolateToBoundaryCellVector(dfloat* phiCell,mesh_t *mesh){

	dfloat* phiFace;
	dfloat* fx;

	createField(&phiFace,surfaceVector,mesh);
	createField(&fx,surfaceScalar,mesh);

	interpolateCellToFaceLinear(phiCell,phiFace,volumeVector,mesh);


	// compute the extrapolation factor

	for (dint face = 0; face < mesh->nBoundaryFaces; face++) {

		dint fid   = mesh->boundaryFaces[face];
		dint cell  =  mesh->faceToOwnerCell[fid];

		dint lface = fid % mesh->ncellFaces; // local face

		dint oface = 0; // opposite face

		if(lface==mesh->leftFace)   {oface=mesh->rightFace;}
		if(lface==mesh->rightFace)  {oface=mesh->leftFace;}
		if(lface==mesh->frontFace)  {oface=mesh->backFace;}
		if(lface==mesh->backFace)   {oface=mesh->frontFace;}
		if(lface==mesh->topFace)    {oface=mesh->bottomFace;}
		if(lface==mesh->bottomFace) {oface=mesh->topFace;}

		dint ofid = oface + cell*mesh->ncellFaces;

		dfloat PN =   mesh->magdf[ofid];
		dfloat bN =   mesh->magdf[fid] + PN;

		fx[fid] = bN/PN;

	}

	// extrapolate to cell face

	for (dint face = 0; face < mesh->nBoundaryFaces; face++) {

		dint fid    = mesh->boundaryFaces[face];
		dint cell   = mesh->faceToOwnerCell[fid];

		dint lface = fid % mesh->ncellFaces; // local face

		dint oface = 0; // opposite face

		if(lface==mesh->leftFace)   {oface=mesh->rightFace;}
		if(lface==mesh->rightFace)  {oface=mesh->leftFace;}
		if(lface==mesh->frontFace)  {oface=mesh->backFace;}
		if(lface==mesh->backFace)   {oface=mesh->frontFace;}
		if(lface==mesh->topFace)    {oface=mesh->bottomFace;}
		if(lface==mesh->bottomFace) {oface=mesh->topFace;}

		dint ofid = oface + cell*mesh->ncellFaces;


		for (int dim = 0; dim < mesh->ndimension; dim++) {

			int fid_dim   = dim + fid*mesh->ndimension;
			int cell_dim  = dim + cell*mesh->ndimension;
			int ofid_dim  = dim + ofid*mesh->ndimension;

			phiFace[fid_dim] = fx[fid]*phiCell[cell_dim] -
					(fx[fid]-1.0)*phiFace[ofid_dim];
		}


	}

	// Copy data to boundary cell

	for (dint face = 0; face < mesh->nBoundaryFaces; face++) {

		dint fid    = mesh->boundaryFaces[face];
		dint cell   = mesh->faceToBoundaryCell[fid];

		for (int dim = 0; dim < mesh->ndimension; dim++) {

			int fid_dim   = dim + fid*mesh->ndimension;
			int cell_dim  = dim + cell*mesh->ndimension;

			phiCell[cell_dim] = phiFace[fid_dim] ;
		}

	}

	free(phiFace);
	free(fx);


//	strcpy(mesh->variableName, "phiCell");
//	printVolumeVector(phiCell, mesh);
}
