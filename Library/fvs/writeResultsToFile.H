#include "mesh.h"

//#include "../Library/fvs/writeResultsToFile.H"

void writeResultsToFile(dfloat* phi, int fieldType, dfloat lengthScaleFactor, dfloat fieldScaleFactor,mesh_t *mesh);

void writeResultsToFile(dfloat* phi, int fieldType, dfloat lengthScaleFactor, dfloat fieldScaleFactor,mesh_t *mesh){


	if (mesh->timeStep % mesh->WriteFrequency != 0) {

		return;
	}

	strcpy(mesh->directoryName, "testCase/");
	mkdir(mesh->directoryName, 0777);

	strcat(mesh->directoryName, "/");

	strcat(mesh->directoryName, mesh->variableName);

	//=================================================================================!

	FILE * fptr;


	strcpy(mesh->fileName, mesh->directoryName); strcat(mesh->fileName,  ".plt");

	if(mesh->timeStep==0){
		fptr = fopen (mesh->fileName, "w");
	} else{
		fptr = fopen (mesh->fileName, "a");
	}
	if(fptr == NULL) {
		printf("Could not create file.\n");
		exit(-1);
	}

	//------------------------------------------------------------------!
	//  Volume Scalar
	//------------------------------------------------------------------!

	if (fieldType == volumeScalar ){

		if(mesh->timeStep==0){
			fprintf(fptr,"VARIABLES = \"X\" \"Z\" \"Y\" \"%s\"\n",mesh->variableName);
		}

		fprintf(fptr,"ZONE T=\"%d\", I=%d, J=%d, K=%d, DATAPACKING=POINT, SOLUTIONTIME=%9.5e,STRANDID=1 \t \n",
				      mesh->timeStep,mesh->ncellsx,mesh->ncellsy,mesh->ncellsz,mesh->time);

		for(dint k = 0; k < mesh->ncellsz; k++){
			for (dint j = 0;j < mesh->ncellsy;j++){
				for (dint i = 0;i < mesh->ncellsx;i++){

					dint cell = computeLinearIndex(i,j,k, mesh->ncellsx, mesh->ncellsy,mesh->ncellsz);

					dint id_x = mesh->dim_x + cell*mesh->ndimension;
					dint id_y = mesh->dim_y + cell*mesh->ndimension;
					dint id_z = mesh->dim_z + cell*mesh->ndimension;

					dfloat X = mesh->C[id_x]*lengthScaleFactor;
					dfloat Y = mesh->C[id_y]*lengthScaleFactor;
					dfloat Z = mesh->C[id_z]*lengthScaleFactor;

					dfloat Phi = phi[cell]*fieldScaleFactor;

					fprintf(fptr,"%9.5e\t%9.5e\t%9.5e\t%9.5e\t", X, Z, Y, Phi);

					fprintf(fptr,"\n");


				}
			}
		}
	}

	//------------------------------------------------------------------!
	//  Volume Vector
	//------------------------------------------------------------------!

	if (fieldType == volumeVector ){

		if(mesh->timeStep==0){
			fprintf(fptr,"VARIABLES = \"X\" \"Z\" \"Y\" \"%sX\" \"%sY\" \"%sZ\"\n",
					      mesh->variableName,mesh->variableName,mesh->variableName);
		}

		fprintf(fptr,"ZONE T=\"%d\", I=%d, J=%d, K=%d, DATAPACKING=POINT, SOLUTIONTIME=%9.5e,STRANDID=1 \t \n",
				      mesh->timeStep,mesh->ncellsx,mesh->ncellsy,mesh->ncellsz,mesh->time);

		for(dint k = 0; k < mesh->ncellsz; k++){
			for (dint j = 0;j < mesh->ncellsy;j++){
				for (dint i = 0;i < mesh->ncellsx;i++){

					dint cell = computeLinearIndex(i,j,k, mesh->ncellsx, mesh->ncellsy,mesh->ncellsz);

					dint id_x = mesh->dim_x + cell*mesh->ndimension;
					dint id_y = mesh->dim_y + cell*mesh->ndimension;
					dint id_z = mesh->dim_z + cell*mesh->ndimension;

					dfloat X = mesh->C[id_x]*lengthScaleFactor;
					dfloat Y = mesh->C[id_y]*lengthScaleFactor;
					dfloat Z = mesh->C[id_z]*lengthScaleFactor;

					dfloat Phi_X = phi[id_x]*fieldScaleFactor;
					dfloat Phi_Y = phi[id_y]*fieldScaleFactor;
					dfloat Phi_Z = phi[id_z]*fieldScaleFactor;

					fprintf(fptr,"%9.5e\t%9.5e\t%9.5e\t%9.5e\t%9.5e\t%9.5e\t", X, Z, Y, Phi_X, Phi_Y, Phi_Z);

					fprintf(fptr,"\n");


				}
			}
		}
	}

	fclose(fptr);


}


