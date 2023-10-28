#include "mesh.h"

//#include "../Library/fvs/meshConstructDestructReconstructFields.H"

void meshConstructReconstructFields(mesh_t *mesh);
void meshDestructReconstructFields(mesh_t *mesh);

void meshConstructReconstructFields(mesh_t *mesh){

	 createField(&mesh->SnDyadSf,surfaceTensor,mesh);
	 createField(&mesh->summedSnDyadSf,volumeTensor,mesh);
	 createField(&mesh->inverseOfSummedSnDyadSf,volumeTensor,mesh);

}

void meshDestructReconstructFields(mesh_t *mesh){

	free(mesh->SnDyadSf);
	free(mesh->summedSnDyadSf);
	free(mesh->inverseOfSummedSnDyadSf);

}
