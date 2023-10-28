#include "mesh.h"
#include "../../setUpCase/FieldNames.H"

//#include "../Library/fvs/constructPEquation.H"

void constructPEquation(mesh_t *mesh);
void setReferencePressureForNullSpace(mesh_t *mesh);


void constructPEquation(mesh_t *mesh){

	reinitializeField(mesh->aP,volumeScalar,mesh);
	reinitializeField(mesh->aN,surfaceScalar,mesh);
	reinitializeField(mesh->bP,volumeScalar,mesh);

	dint Peqn;

	Peqn = fvmLaplacianOfScalarWithVectorCoefficient (1.0, mesh->rAPf, mesh->p,p,mesh);
	Peqn = fvcDivergenceOfFlux(-1.0,mesh->phiHbyA,mesh);




//	setReferencePressureForNullSpace(mesh);



}



void setReferencePressureForNullSpace(mesh_t *mesh){

	int cell = 0; // Reference cell

	mesh->aP[cell] = 1; // Make the diagonal value to be 1
	mesh->bP[cell] = 1; // Reference pressure value


	for (int face = 0; face < mesh->ncellFaces; face++) {

			int fid 	= face+cell*mesh->ncellFaces;
			int cell_n 	= mesh->faceToNeighborCell[fid];

			mesh->aN[fid] = 0;  // Make the off diagonal value to be 0

	}


}
