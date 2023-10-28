/******************************************************************************
   Grid Generator
*******************************************************************************/

//valgrind --tool=memcheck --leak-check=yes --show-reachable=yes --num-callers=20 --track-fds=yes ./main

#include "mesh.h"

int main(int argc, char **argv){

	clock_t t;

	t = clock();

	mesh_t *mesh = (mesh_t*) calloc(1, sizeof(mesh_t));

	setMeshControls(mesh);

	assignMeshDimensions(mesh);

	calculateIndices(mesh);

	createField(&mesh->vertices,vertexVector,mesh);

	printf("Generating Boundary Mesh  \n");

	setMeshCorners(mesh);

	assignMeshCorners(mesh);

	meshConstructComputationalDomainFields(mesh);

	generateComputationalDomain(mesh);

	generateEdgeMesh(mesh);

	printf("Generating Surface Mesh \n");

	createSurfaceMesh(mesh);

	createVolumeMesh(mesh);

	convertCylindricalPolarToCartesian(mesh);

	convertSphericalPolarToCartesian(mesh);

	printf("Writing Surface Mesh \n");

	writeFiniteDifferenceMesh(mesh);

	printf("Writing Finite Volume Surface Mesh \n");

	writeFiniteVolumeMesh(mesh);

	meshDestructComputationalDomainFields(mesh);

	free(mesh->vertices);

	free(mesh);

	t = clock() - t;
	double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds

	printf("Mesh Generation took %f seconds to execute \n", time_taken);

}
