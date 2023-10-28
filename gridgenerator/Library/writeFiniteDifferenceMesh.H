#include "mesh.h"

void writeFiniteDifferenceMesh(mesh_t *mesh);

void writeFiniteDifferenceMesh(mesh_t *mesh){


	int n = 50;
    char *fileName ;
    fileName	= (char*) calloc(n, sizeof(char));
    strcpy(fileName, "finiteDifferenceMesh");

    // create directory
	char *directory ;
	directory		= (char*) calloc(n, sizeof(char));
    strcpy(directory, "Mesh/");
    mkdir(directory, 0777);


	strcat(directory, fileName);


    char *newFilelName ;
    newFilelName	= (char*) calloc(n, sizeof(char));

    FILE * fptr;

    strcpy(newFilelName, directory); strcat(newFilelName, ".plt");


    fptr = fopen (newFilelName, "w");
    if(fptr == NULL) {
    	printf("Could not create file.\n");
    	exit(-1);
    }

	fprintf(fptr,"VARIABLES = \"X\" \"Y\" \"Z\" \n");

	fprintf(fptr,"ZONE I=%d, J=%d, K=%d, DATAPACKING=POINT \t \n",mesh->nverticesx,mesh->nverticesy,mesh->nverticesz);


	for(dint k = 0; k < mesh->nverticesz; k++){
		for (dint j = 0;j < mesh->nverticesy;j++){
			for (dint i = 0;i < mesh->nverticesx;i++){

				dint vertex = computeLinearIndex(i,j,k, mesh->nverticesx, mesh->nverticesy, mesh->nverticesz);

				dint vid_x = mesh->dim_x + vertex*mesh->ndimension;
				dint vid_y = mesh->dim_y + vertex*mesh->ndimension;
				dint vid_z = mesh->dim_z + vertex*mesh->ndimension;

				fprintf(fptr,"%lf\t%lf\t%lf\t",mesh->vertices[vid_x],mesh->vertices[vid_y],mesh->vertices[vid_z]);

				fprintf(fptr,"\n");


			}
		}
	}


    fclose(fptr);



    free(newFilelName);
	free(directory);
    free(fileName);



}
