#include "mesh.h"

//#include "../setUpCase/setReferenceHeight.H"

void setReferenceHeight(mesh_t *mesh);

void setReferenceHeight(mesh_t *mesh){

	mesh->Xref[mesh->dim_x] = 0.00;      // x value
	mesh->Xref[mesh->dim_y] = 0.00;      // y value
	mesh->Xref[mesh->dim_z] = 1.00;      // z value

}

