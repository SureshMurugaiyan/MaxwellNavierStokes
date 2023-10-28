/******************************************************************************
   Maxwell-Navier-Stokes Solver Fusion Liquid Wall
*******************************************************************************/

//valgrind --tool=memcheck --leak-check=yes --show-reachable=yes --num-callers=20 --track-fds=yes ./main
// clear && make PETSC_DIR="/home/sureshm/Softwares/petsc/"  PETSC_ARCH="arch-linux-c-debug" && ./main

#include "mesh.h"
#include "petscksp.h"

static char help[] = "Solves Maxwell-Navier-Stokes Equation.\n\n";

int main(int argc, char **argv){

	clock_t time;

	time = clock();

	mesh_t *mesh = (mesh_t*) calloc(1, sizeof(mesh_t));

	mesh->CPUstarttime = time;

	PetscErrorCode ierr;
	PetscMPIInt    size;
	PetscInt 	   n ;

	ierr = PetscInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
	if (size != 1) SETERRQ(PETSC_COMM_WORLD,1,"This is a uniprocessor Solver!");
	ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);

	meshReadStatistics(mesh);

	meshSettings(mesh);

	meshConstructVerticesFields(mesh);

	meshReadVertices(mesh);

	meshConstructConnectivityFields(mesh);

	meshGenerateCellConnectivity(mesh);

	meshCalculateNoOfBoundaryandInnerFaces(mesh);

	meshConstructFaceToCellConnectivityFields(mesh);

	meshGenerateFaceToCellConnectivity(mesh);

	meshConstructMetricsFields(mesh);

	meshGenerateCellMetrics(mesh);

	meshGenerateBoundaryPatchID(mesh);

	meshCountNoOfBoundaryPatches(mesh);

	meshConstructPatchConnectivityFields(mesh);

	meshGenerateBoundaryPatchConnectivity(mesh);

	meshSetUpPatchType(mesh);

	meshSetUpPeriodicNeighborPatches(mesh);

	meshGeneratePeriodicConnectivity(mesh);

	meshConstructInnerAndPeriodicFacesFields(mesh);

	meshCalculateInnerFaceIDIncludingPeriodicFaces(mesh);

	meshConstructDomainFields(mesh);

	setUpDomains(mesh);

	printf("\nWriting Mesh and Mesh Metrics \n");

	writeMesh(mesh);

	meshConstructReconstructFields(mesh);

	generateMetricsForFluxReconstruction(mesh);

	mainNavierStokesSolve(mesh);






	meshDestructReconstructFields(mesh);


	meshDestructDomainFields(mesh);

	meshDestructInnerAndPeriodicFacesFields(mesh);

	meshDestructPatchConnectivityFields(mesh);

	meshDestructMetricsFields(mesh);

	meshDestructFaceToCellConnectivityFields(mesh);

	meshDestructConnectivityFields(mesh);

	meshDestructVerticesFields(mesh);

	PetscCall(PetscFinalize());

	time = clock() - time;

	dfloat time_taken = ((dfloat)time)/CLOCKS_PER_SEC; // in seconds

	printf("Solver took %f seconds to execute \n", time_taken);

	free(mesh);

}
