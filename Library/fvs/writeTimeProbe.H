#include "mesh.h"

//#include "../Library/fvs/writeTimeProbe.H"

void writeTimeProbe(dfloat *phi,int fieldType,mesh_t *mesh);

void writeTimeProbe(dfloat *phi,int fieldType,mesh_t *mesh){

	if (mesh->timeStep % mesh->WriteFrequency != 0) {

		return;
	}

    strcpy(mesh->directoryName, "testCaseTimeProbe/");
    mkdir(mesh->directoryName, 0777);

	strcat(mesh->directoryName, mesh->variableName);

	dint ndimension   = getNdimension(fieldType);

	//=================================================================================!

	FILE * fptr;

	for( int dim = 0; dim < ndimension;dim++) {

		if ((dim == 0) && (ndimension==1)) { strcpy(mesh->fileName, mesh->directoryName); strcat(mesh->fileName,  ".txt");};
		if ((dim == 0) && (ndimension==3)) { strcpy(mesh->fileName, mesh->directoryName); strcat(mesh->fileName, "X.txt");};
		if ((dim == 1) && (ndimension==3)) { strcpy(mesh->fileName, mesh->directoryName); strcat(mesh->fileName, "Y.txt");};
		if ((dim == 2) && (ndimension==3)) { strcpy(mesh->fileName, mesh->directoryName); strcat(mesh->fileName, "Z.txt");};



		if(mesh->timeStep==0){

			fptr = fopen (mesh->fileName, "w");
		}else{

			fptr = fopen (mesh->fileName, "a");
		}

		if(fptr == NULL) {
			printf("Could not create file.\n");
			exit(-1);
		}


		fprintf(fptr,"%d \t",mesh->timeStep);
		fprintf(fptr,"%9.10e \t",mesh->time);
		fprintf(fptr,"%9.10e",phi[dim]);
		fprintf(fptr,"\n");



		fclose(fptr);
	}



}
