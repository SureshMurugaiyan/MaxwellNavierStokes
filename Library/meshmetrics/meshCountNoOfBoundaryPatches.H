#include "mesh.h"

void meshCountNoOfBoundaryPatches(mesh_t *mesh);

void meshCountNoOfBoundaryPatches(mesh_t *mesh){

	int count = 0;
	int j = 0;

	// Use nested loop
	// Pick elements one by one
	for (int i = 0; i < mesh->nBoundaryFaces; i++) {
		// check if the element has appeared to the left side of array
		for (j = 0; j < i; j++){
			// if the elements has already appeared, break the for loop
			// and pick the next element
			if (mesh->boundaryFacePatchID[i] == mesh->boundaryFacePatchID[j])
			{
				break;
			}
		}

		// if the elements has already not appeared
		// then you will reach the condition i=j
		// increase the count by one
		if (i == j) {count++;}
	}

	mesh->npatchs  = count;



//	printf ( "!--------------------------------------------- !\n" );
//	printf ( "! Number of Boundary Patches                   !\n" );
//	printf ( "!--------------------------------------------- !\n!" );
//
//	printf ( " %5d \n!", mesh->npatchs  );
//	printf ( "--------------------------------------------- !\n" );

}
