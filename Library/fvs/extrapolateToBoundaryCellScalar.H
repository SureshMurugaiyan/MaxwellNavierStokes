#include "mesh.h"

//#include "../Library/fvs/extrapolateToBoundaryCellScalar.H"


void extrapolateToBoundaryCellScalar(dfloat* phiCell,mesh_t *mesh);

void extrapolateToBoundaryCellScalar(dfloat* phiCell,mesh_t *mesh){

	dfloat* phiFace;
	dfloat* fx;

	createField(&phiFace,surfaceScalar,mesh);
	createField(&fx,surfaceScalar,mesh);

	interpolateCellToFaceLinear(phiCell,phiFace,volumeScalar,mesh);

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


		phiFace[fid] = fx[fid]*phiCell[cell] -(fx[fid]-1.0)*phiFace[ofid];


	}

	// Copy data to boundary cell

	for (dint face = 0; face < mesh->nBoundaryFaces; face++) {

		dint fid    = mesh->boundaryFaces[face];
		dint cell   = mesh->faceToBoundaryCell[fid];

		phiCell[cell] = phiFace[fid] ;

	}


	free(phiFace);
	free(fx);


//	strcpy(mesh->variableName, "phiCell");
//	printVolumeScalar(phiCell, mesh);
}
