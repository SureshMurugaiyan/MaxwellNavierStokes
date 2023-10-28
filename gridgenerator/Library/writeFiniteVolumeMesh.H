#include "mesh.h"

void writeFiniteVolumeMesh(mesh_t *mesh);

void writeFiniteVolumeMesh(mesh_t *mesh){



	mesh->cellToVertex 			= (dint*) calloc(mesh->ncells*mesh->ncellVertices, sizeof(dint));
	for (int i = 0; i < (mesh->ncells*mesh->ncellVertices); i++) {
		mesh->cellToVertex[i] = 0;
	}


	meshCellToVerticesConnectivity(mesh);

    // create directory
	int n = 50;

	char *directory ;
	directory		= (char*) calloc(n, sizeof(char));
    strcpy(directory, "Mesh/");
    mkdir(directory, 0777);
    free(directory);


    // create fileName
	char *fileName ;
	fileName	= (char*) calloc(n, sizeof(char));
	strcpy(fileName, "Mesh/finiteVolumeMesh.txt");


	FILE * fptr;

	fptr = fopen (fileName, "w");
	if(fptr == NULL) {
		printf("Could not create file.\n");
		exit(-1);
	}

	for( int cell = 0; cell < mesh->ncells;cell++){
		for(int vertex = 0;vertex < mesh->ncellVertices;vertex++){
			for(int dimension = 0;dimension < mesh->ndimension;dimension++){

				int id 	= vertex + cell*mesh->ncellVertices ;
				int vid = mesh->cellToVertex[id];
				int vid_dim = dimension + vid*mesh->ndimension ;

				fprintf(fptr,"%lf \n",mesh->vertices[vid_dim]);

			}
		}
	}
	fclose(fptr);



    // Write Mesh Statistics
    strcpy(fileName, "Mesh/finiteVolumeMeshStatistics.txt");
    fptr = fopen (fileName, "w");
    if(fptr == NULL) {
    	printf("Could not create file.\n");
    	exit(-1);
    }

    fprintf(fptr,"%d \n",mesh->ncellsx);
    fprintf(fptr,"%d \n",mesh->ncellsy);
    fprintf(fptr,"%d \n",mesh->ncellsz);

    fclose(fptr);




    free(fileName);
	free(mesh->cellToVertex);

}
