#include "mesh.h"
#include "BoundaryNames.H"
#include "FieldNames.H"

//#include "../setUpCase/getBoundaryConditionValue.H"

void getBoundaryConditionValue(dint patch, dint phiName,  dfloat* bcValue, mesh_t *mesh);

void getBoundaryConditionValue(dint patch, dint phiName,  dfloat* bcValue, mesh_t *mesh){

	dfloat phiVal = 0.01;

	dint nPatchFaces = mesh->npatchFaces[patch];
	dint face = mesh->patchStartFaceIndex[patch];

	dint patchID = mesh->patchID[patch];

	dint dimx = 0; // dim = 0 gives x component
	dint dimy = 1; // dim = 1 gives y component
	dint dimz = 2; // dim = 2 gives z component

	dint PatchFacex,PatchFacey,PatchFacez;
	dint fidx,fidy,fidz;
	dfloat x,y,z;


	switch (patchID)
	{
//====================================================================================================!
	case topBoundary:
		switch (phiName)
		{
		case U:
			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++){

				dint fid  = mesh->patchToBoundaryFace[face];

				fidx	= dimx+fid*mesh->ndimension;
				fidy	= dimy+fid*mesh->ndimension;
				fidz	= dimz+fid*mesh->ndimension;

				x = mesh->Cf[fidx];
				y = mesh->Cf[fidy];
				z = mesh->Cf[fidz];


				dint idx	= dimx+PatchFace*mesh->ndimension;
				dint idy	= dimy+PatchFace*mesh->ndimension;
				dint idz	= dimz+PatchFace*mesh->ndimension;


				bcValue[idx]	= 1.0;
				bcValue[idy]	= 0.0;
				bcValue[idz]	= 0.0;

				face++;
			}
			break;
		case p:
			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++){

				dint fid  = mesh->patchToBoundaryFace[face];

				fidx	= dimx+fid*mesh->ndimension;
				fidy	= dimy+fid*mesh->ndimension;
				fidz	= dimz+fid*mesh->ndimension;

				x = mesh->Cf[fidx];
				y = mesh->Cf[fidy];
				z = mesh->Cf[fidz];

				bcValue[PatchFace]	= 0.0;

				face++;
			}
			break;
		case alpha:
			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++){

				dint fid  = mesh->patchToBoundaryFace[face];

				fidx	= dimx+fid*mesh->ndimension;
				fidy	= dimy+fid*mesh->ndimension;
				fidz	= dimz+fid*mesh->ndimension;

				x = mesh->Cf[fidx];
				y = mesh->Cf[fidy];
				z = mesh->Cf[fidz];


				bcValue[PatchFace]	= 0.0;

				face++;
			}
			break;
		case phi:
			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++){

				dint fid  = mesh->patchToBoundaryFace[face];

				fidx	= dimx+fid*mesh->ndimension;
				fidy	= dimy+fid*mesh->ndimension;
				fidz	= dimz+fid*mesh->ndimension;

				x = mesh->Cf[fidx];
				y = mesh->Cf[fidy];
				z = mesh->Cf[fidz];


				bcValue[PatchFace]	= 0.0;

				face++;
			}
			break;
		case A:
			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++){

				dint fid  = mesh->patchToBoundaryFace[face];

				fidx	= dimx+fid*mesh->ndimension;
				fidy	= dimy+fid*mesh->ndimension;
				fidz	= dimz+fid*mesh->ndimension;

				x = mesh->Cf[fidx];
				y = mesh->Cf[fidy];
				z = mesh->Cf[fidz];


				dint idx	= dimx+PatchFace*mesh->ndimension;
				dint idy	= dimy+PatchFace*mesh->ndimension;
				dint idz	= dimz+PatchFace*mesh->ndimension;


				bcValue[idx]	= 0.0;
				bcValue[idy]	= 0.0;
				bcValue[idz]	= 0.0;

				face++;
			}
			break;

		}
		break;
//====================================================================================================!
	case bottomBoundary:
		switch (phiName)
		{
		case U:
			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++){

				dint fid  = mesh->patchToBoundaryFace[face];

				fidx	= dimx+fid*mesh->ndimension;
				fidy	= dimy+fid*mesh->ndimension;
				fidz	= dimz+fid*mesh->ndimension;

				x = mesh->Cf[fidx];
				y = mesh->Cf[fidy];
				z = mesh->Cf[fidz];


				dint idx	= dimx+PatchFace*mesh->ndimension;
				dint idy	= dimy+PatchFace*mesh->ndimension;
				dint idz	= dimz+PatchFace*mesh->ndimension;


				bcValue[idx]	= 0.0;
				bcValue[idy]	= 0.0;
				bcValue[idz]	= 0.0;

				face++;
			}
			break;
		case p:
			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++){

				dint fid  = mesh->patchToBoundaryFace[face];

				fidx	= dimx+fid*mesh->ndimension;
				fidy	= dimy+fid*mesh->ndimension;
				fidz	= dimz+fid*mesh->ndimension;

				x = mesh->Cf[fidx];
				y = mesh->Cf[fidy];
				z = mesh->Cf[fidz];

				bcValue[PatchFace]	= 0.0;

				face++;
			}
			break;
		case alpha:
			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++){

				dint fid  = mesh->patchToBoundaryFace[face];

				fidx	= dimx+fid*mesh->ndimension;
				fidy	= dimy+fid*mesh->ndimension;
				fidz	= dimz+fid*mesh->ndimension;

				x = mesh->Cf[fidx];
				y = mesh->Cf[fidy];
				z = mesh->Cf[fidz];


				bcValue[PatchFace]	= 0.0;

				face++;
			}
			break;
		case phi:
			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++){

				dint fid  = mesh->patchToBoundaryFace[face];

				fidx	= dimx+fid*mesh->ndimension;
				fidy	= dimy+fid*mesh->ndimension;
				fidz	= dimz+fid*mesh->ndimension;

				x = mesh->Cf[fidx];
				y = mesh->Cf[fidy];
				z = mesh->Cf[fidz];


				bcValue[PatchFace]	= 0.0;

				face++;
			}
			break;
		case A:
			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++){

				dint fid  = mesh->patchToBoundaryFace[face];

				fidx	= dimx+fid*mesh->ndimension;
				fidy	= dimy+fid*mesh->ndimension;
				fidz	= dimz+fid*mesh->ndimension;

				x = mesh->Cf[fidx];
				y = mesh->Cf[fidy];
				z = mesh->Cf[fidz];


				dint idx	= dimx+PatchFace*mesh->ndimension;
				dint idy	= dimy+PatchFace*mesh->ndimension;
				dint idz	= dimz+PatchFace*mesh->ndimension;


				bcValue[idx]	= 0.0;
				bcValue[idy]	= 0.0;
				bcValue[idz]	= 0.0;

				face++;
			}
			break;

		}
		break;
//====================================================================================================!
	case frontBoundary:
		switch (phiName)
		{
		case U:
			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++){

				dint fid  = mesh->patchToBoundaryFace[face];

				fidx	= dimx+fid*mesh->ndimension;
				fidy	= dimy+fid*mesh->ndimension;
				fidz	= dimz+fid*mesh->ndimension;

				x = mesh->Cf[fidx];
				y = mesh->Cf[fidy];
				z = mesh->Cf[fidz];


				dint idx	= dimx+PatchFace*mesh->ndimension;
				dint idy	= dimy+PatchFace*mesh->ndimension;
				dint idz	= dimz+PatchFace*mesh->ndimension;


				bcValue[idx]	= 0.0;
				bcValue[idy]	= 0.0;
				bcValue[idz]	= 0.0;

				face++;
			}
			break;
		case p:
			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++){

				dint fid  = mesh->patchToBoundaryFace[face];

				fidx	= dimx+fid*mesh->ndimension;
				fidy	= dimy+fid*mesh->ndimension;
				fidz	= dimz+fid*mesh->ndimension;

				x = mesh->Cf[fidx];
				y = mesh->Cf[fidy];
				z = mesh->Cf[fidz];

				bcValue[PatchFace]	= 0.0;

				face++;
			}
			break;
		case alpha:
			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++){

				dint fid  = mesh->patchToBoundaryFace[face];

				fidx	= dimx+fid*mesh->ndimension;
				fidy	= dimy+fid*mesh->ndimension;
				fidz	= dimz+fid*mesh->ndimension;

				x = mesh->Cf[fidx];
				y = mesh->Cf[fidy];
				z = mesh->Cf[fidz];


				bcValue[PatchFace]	= 0.0;

				face++;
			}
			break;
		case phi:
			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++){

				dint fid  = mesh->patchToBoundaryFace[face];

				fidx	= dimx+fid*mesh->ndimension;
				fidy	= dimy+fid*mesh->ndimension;
				fidz	= dimz+fid*mesh->ndimension;

				x = mesh->Cf[fidx];
				y = mesh->Cf[fidy];
				z = mesh->Cf[fidz];


				bcValue[PatchFace]	=  phiVal;

				face++;
			}
			break;
		case A:
			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++){

				dint fid  = mesh->patchToBoundaryFace[face];

				fidx	= dimx+fid*mesh->ndimension;
				fidy	= dimy+fid*mesh->ndimension;
				fidz	= dimz+fid*mesh->ndimension;

				x = mesh->Cf[fidx];
				y = mesh->Cf[fidy];
				z = mesh->Cf[fidz];


				dint idx	= dimx+PatchFace*mesh->ndimension;
				dint idy	= dimy+PatchFace*mesh->ndimension;
				dint idz	= dimz+PatchFace*mesh->ndimension;


				bcValue[idx]	= 0.0;
				bcValue[idy]	= 0.0;
				bcValue[idz]	= 0.0;

				face++;
			}
			break;


		}
		break;
//====================================================================================================!
	case backBoundary:
		switch (phiName)
		{
		case U:
			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++){

				dint fid  = mesh->patchToBoundaryFace[face];

				fidx	= dimx+fid*mesh->ndimension;
				fidy	= dimy+fid*mesh->ndimension;
				fidz	= dimz+fid*mesh->ndimension;

				x = mesh->Cf[fidx];
				y = mesh->Cf[fidy];
				z = mesh->Cf[fidz];


				dint idx	= dimx+PatchFace*mesh->ndimension;
				dint idy	= dimy+PatchFace*mesh->ndimension;
				dint idz	= dimz+PatchFace*mesh->ndimension;


				bcValue[idx]	= 0.0;
				bcValue[idy]	= 0.0;
				bcValue[idz]	= 0.0;

				face++;
			}
			break;
		case p:
			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++){

				dint fid  = mesh->patchToBoundaryFace[face];

				fidx	= dimx+fid*mesh->ndimension;
				fidy	= dimy+fid*mesh->ndimension;
				fidz	= dimz+fid*mesh->ndimension;

				x = mesh->Cf[fidx];
				y = mesh->Cf[fidy];
				z = mesh->Cf[fidz];

				bcValue[PatchFace]	= 0.0;

				face++;
			}
			break;
		case alpha:
			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++){

				dint fid  = mesh->patchToBoundaryFace[face];

				fidx	= dimx+fid*mesh->ndimension;
				fidy	= dimy+fid*mesh->ndimension;
				fidz	= dimz+fid*mesh->ndimension;

				x = mesh->Cf[fidx];
				y = mesh->Cf[fidy];
				z = mesh->Cf[fidz];


				bcValue[PatchFace]	= 0.0;

				face++;
			}
			break;
		case phi:
			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++){

				dint fid  = mesh->patchToBoundaryFace[face];

				fidx	= dimx+fid*mesh->ndimension;
				fidy	= dimy+fid*mesh->ndimension;
				fidz	= dimz+fid*mesh->ndimension;

				x = mesh->Cf[fidx];
				y = mesh->Cf[fidy];
				z = mesh->Cf[fidz];


				bcValue[PatchFace]	= -phiVal;

				face++;
			}
			break;
		case A:
			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++){

				dint fid  = mesh->patchToBoundaryFace[face];

				fidx	= dimx+fid*mesh->ndimension;
				fidy	= dimy+fid*mesh->ndimension;
				fidz	= dimz+fid*mesh->ndimension;

				x = mesh->Cf[fidx];
				y = mesh->Cf[fidy];
				z = mesh->Cf[fidz];


				dint idx	= dimx+PatchFace*mesh->ndimension;
				dint idy	= dimy+PatchFace*mesh->ndimension;
				dint idz	= dimz+PatchFace*mesh->ndimension;


				bcValue[idx]	= 0.0;
				bcValue[idy]	= 0.0;
				bcValue[idz]	= 0.0;

				face++;
			}
			break;
		}
		break;

//====================================================================================================!














//====================================================================================================!
	case leftBoundary:
		switch (phiName)
		{
		case U:
			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++){

				dint fid  = mesh->patchToBoundaryFace[face];

				fidx	= dimx+fid*mesh->ndimension;
				fidy	= dimy+fid*mesh->ndimension;
				fidz	= dimz+fid*mesh->ndimension;

				x = mesh->Cf[fidx];
				y = mesh->Cf[fidy];
				z = mesh->Cf[fidz];


				dint idx	= dimx+PatchFace*mesh->ndimension;
				dint idy	= dimy+PatchFace*mesh->ndimension;
				dint idz	= dimz+PatchFace*mesh->ndimension;


				bcValue[idx]	= 0.0;
				bcValue[idy]	= 0.0;
				bcValue[idz]	= 0.0;

				face++;
			}
			break;
		case p:
			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++){

				dint fid  = mesh->patchToBoundaryFace[face];

				fidx	= dimx+fid*mesh->ndimension;
				fidy	= dimy+fid*mesh->ndimension;
				fidz	= dimz+fid*mesh->ndimension;

				x = mesh->Cf[fidx];
				y = mesh->Cf[fidy];
				z = mesh->Cf[fidz];

				bcValue[PatchFace]	= 0.0;

				face++;
			}
			break;
		case alpha:
			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++){

				dint fid  = mesh->patchToBoundaryFace[face];

				fidx	= dimx+fid*mesh->ndimension;
				fidy	= dimy+fid*mesh->ndimension;
				fidz	= dimz+fid*mesh->ndimension;

				x = mesh->Cf[fidx];
				y = mesh->Cf[fidy];
				z = mesh->Cf[fidz];


				bcValue[PatchFace]	= 0.0;

				face++;
			}
			break;
		case phi:
			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++){

				dint fid  = mesh->patchToBoundaryFace[face];

				fidx	= dimx+fid*mesh->ndimension;
				fidy	= dimy+fid*mesh->ndimension;
				fidz	= dimz+fid*mesh->ndimension;

				x = mesh->Cf[fidx];
				y = mesh->Cf[fidy];
				z = mesh->Cf[fidz];


				bcValue[PatchFace]	= 0.0;

				face++;
			}
			break;
		case A:
			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++){

				dint fid  = mesh->patchToBoundaryFace[face];

				fidx	= dimx+fid*mesh->ndimension;
				fidy	= dimy+fid*mesh->ndimension;
				fidz	= dimz+fid*mesh->ndimension;

				x = mesh->Cf[fidx];
				y = mesh->Cf[fidy];
				z = mesh->Cf[fidz];


				dint idx	= dimx+PatchFace*mesh->ndimension;
				dint idy	= dimy+PatchFace*mesh->ndimension;
				dint idz	= dimz+PatchFace*mesh->ndimension;


				bcValue[idx]	= 0.0;
				bcValue[idy]	= 0.0;
				bcValue[idz]	= 0.0;

				face++;
			}
			break;



		}
		break;
//====================================================================================================!
	case rightBoundary:
		switch (phiName)
		{
		case U:
			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++){

				dint fid  = mesh->patchToBoundaryFace[face];

				fidx	= dimx+fid*mesh->ndimension;
				fidy	= dimy+fid*mesh->ndimension;
				fidz	= dimz+fid*mesh->ndimension;

				x = mesh->Cf[fidx];
				y = mesh->Cf[fidy];
				z = mesh->Cf[fidz];


				dint idx	= dimx+PatchFace*mesh->ndimension;
				dint idy	= dimy+PatchFace*mesh->ndimension;
				dint idz	= dimz+PatchFace*mesh->ndimension;


				bcValue[idx]	= 0.0;
				bcValue[idy]	= 0.0;
				bcValue[idz]	= 0.0;

				face++;
			}
			break;
		case p:
			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++){

				dint fid  = mesh->patchToBoundaryFace[face];

				fidx	= dimx+fid*mesh->ndimension;
				fidy	= dimy+fid*mesh->ndimension;
				fidz	= dimz+fid*mesh->ndimension;

				x = mesh->Cf[fidx];
				y = mesh->Cf[fidy];
				z = mesh->Cf[fidz];

				bcValue[PatchFace]	= 0.0;

				face++;
			}
			break;
		case alpha:
			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++){

				dint fid  = mesh->patchToBoundaryFace[face];

				fidx	= dimx+fid*mesh->ndimension;
				fidy	= dimy+fid*mesh->ndimension;
				fidz	= dimz+fid*mesh->ndimension;

				x = mesh->Cf[fidx];
				y = mesh->Cf[fidy];
				z = mesh->Cf[fidz];


				bcValue[PatchFace]	= 0.0;

				face++;
			}
			break;
		case phi:
			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++){

				dint fid  = mesh->patchToBoundaryFace[face];

				fidx	= dimx+fid*mesh->ndimension;
				fidy	= dimy+fid*mesh->ndimension;
				fidz	= dimz+fid*mesh->ndimension;

				x = mesh->Cf[fidx];
				y = mesh->Cf[fidy];
				z = mesh->Cf[fidz];


				bcValue[PatchFace]	= 0.0;

				face++;
			}
			break;
		case A:
			for (dint PatchFace = 0; PatchFace < nPatchFaces; PatchFace++){

				dint fid  = mesh->patchToBoundaryFace[face];

				fidx	= dimx+fid*mesh->ndimension;
				fidy	= dimy+fid*mesh->ndimension;
				fidz	= dimz+fid*mesh->ndimension;

				x = mesh->Cf[fidx];
				y = mesh->Cf[fidy];
				z = mesh->Cf[fidz];


				dint idx	= dimx+PatchFace*mesh->ndimension;
				dint idy	= dimy+PatchFace*mesh->ndimension;
				dint idz	= dimz+PatchFace*mesh->ndimension;


				bcValue[idx]	= 0.0;
				bcValue[idy]	= 0.0;
				bcValue[idz]	= 0.0;

				face++;
			}
			break;
		}
		break;
	}



}
