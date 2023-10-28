
#include "mesh.h"

void meshReadVertices(mesh_t *mesh);

void meshReadVertices(mesh_t *mesh){


	FILE *myFile;
	myFile = fopen("setUpCase/InputMesh/finiteVolumeMesh.txt", "r");

	if (myFile == NULL){
		printf("Error Reading File 'setUpCase/InputMesh/finiteVolumeMesh.txt' \n");
		printf("File Not Found \n");
		exit (0);
	}


	// Load the inner mesh vertices from the Mesh File

	for (dint i = 0; i < mesh->nivertices*mesh->ndimension; i++){

    	if (fscanf(myFile, "%lf",&mesh->vertices[i] ) > 0){

    	}else{
    		printf("cannot read data:Index exceeded:setUpCase/finiteVolumeMesh.txt");
    	}

	}


	fclose(myFile);



//	strcpy(mesh->variableName, "V");
//	printVertexVector(mesh->vertices, mesh);
}
