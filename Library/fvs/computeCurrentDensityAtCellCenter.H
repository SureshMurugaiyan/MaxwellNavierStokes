#include "mesh.h"

//#include "../Library/fvs/computeCurrentDensityAtCellCenter.H"

void computeCurrentDensityAtCellCenter(mesh_t *mesh);

void computeCurrentDensityAtCellCenter(mesh_t *mesh){

	reconstructFluxFromCellFaceToCellCenter(mesh->sigmagradPhiFlux,mesh->Jphi,mesh);
	reconstructFluxFromCellFaceToCellCenter(mesh->phiJ,mesh->J,mesh);


	// Computing current densities at boundaries
	// Ju at boundaries is already known
	// Need to calculate Jphi at boundaries


	for (dint face = 0; face < mesh->nBoundaryFaces; face++) {

		dint fid    = mesh->boundaryFaces[face];
		dint Bcell  =  mesh->faceToBoundaryCell[fid];
		dint Ocell  =  mesh->faceToOwnerCell[fid];

		for (int dim = 0; dim < mesh->ndimension; dim++) {

			int fid_dim    = dim + fid*mesh->ndimension;
			int Bcell_dim  = dim + Bcell*mesh->ndimension;
			int Ocell_dim  = dim + Ocell*mesh->ndimension;

			mesh->J[Bcell_dim] = - mesh->sigma[Bcell]*(mesh->phi[Bcell] - mesh->phi[Ocell])*mesh->rdelta[fid]*mesh->Sn[fid_dim]
								 + mesh->Ju[Bcell_dim];

		}
	}



//	dfloat *phif;     // electric potential faceValues
//	dfloat *gradPhi;  // gradient of electric potential
//
//
//	createField(&phif,surfaceScalar,mesh);
//	createField(&gradPhi,volumeVector,mesh);
//
//
//	interpolateCellToFaceLinear(mesh->phi,phif,volumeScalar,mesh);
//	fvdGradientOfScalar(-1.0,phif,mesh,gradPhi);
//
//	for (int cell = 0; cell < mesh->nicells; cell++) {
//		for (int dim = 0; dim < mesh->ndimension; dim++) {
//
//			dint cell_dim = dim + cell*mesh->ndimension;
//
//			mesh->Jphi[cell_dim] = mesh->sigma[cell]*gradPhi[cell_dim];
//		}
//	}
//
//
//
//	for (int cell = 0; cell < mesh->nicells; cell++) {
//		for (int dim = 0; dim < mesh->ndimension; dim++) {
//
//			dint cell_dim = dim + cell*mesh->ndimension;
//
//			mesh->J[cell_dim] = mesh->Jphi[cell_dim] +  mesh->Ju[cell_dim] ;
//		}
//	}
//
//	free(phif);
//	free(gradPhi);



}

