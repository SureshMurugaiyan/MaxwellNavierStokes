#include "mesh.h"

//#include "../Library/fvs/writeZLineDataToFile.H"

void writeZLineDataToFile(dfloat* Phi, int fieldType, dint I, dfloat lengthScaleFactor, dfloat fieldScaleFactor,mesh_t *mesh);

void writeZLineDataToFile(dfloat* phi, int fieldType, dint I, dfloat lengthScaleFactor, dfloat fieldScaleFactor,mesh_t *mesh){

	if (I > (mesh->ncellsx-1)) {

		return;
	}


	if (mesh->timeStep % mesh->WriteFrequency != 0) {

		return;
	}



	strcpy(mesh->directoryName, "testCase/");
	mkdir(mesh->directoryName, 0777);

    sprintf(mesh->subdirectoryName,"%d",mesh->timeStep);
    strcat(mesh->directoryName, mesh->subdirectoryName);
	mkdir(mesh->directoryName, 0777);

	strcat(mesh->directoryName, "/");

	strcat(mesh->directoryName, mesh->variableName);

	//=================================================================================!

	FILE * fptr;


	strcpy(mesh->fileName, mesh->directoryName); strcat(mesh->fileName,  "_ZLINE.txt");

	fptr = fopen (mesh->fileName, "w");
	if(fptr == NULL) {
		printf("Could not create file.\n");
		exit(-1);
	}

	dint J = 0; // only one cell in Phi direction

	//------------------------------------------------------------------!
	//  Volume Scalar
	//------------------------------------------------------------------!

	if (fieldType == volumeScalar ){

		fprintf(fptr,"VARIABLES = \"X\" \"Y\" \"Z\" \"%s\"\n",mesh->variableName);

		fprintf(fptr,"ZONE I=%d, DATAPACKING=POINT \t \n",mesh->ncellsx);

		dint j = J;
		dint i = I;

		for (dint k = 0;k < mesh->ncellsz;k++){

			dint cell = computeLinearIndex(i,j,k, mesh->ncellsx, mesh->ncellsy,mesh->ncellsz);

			dint id_x = mesh->dim_x + cell*mesh->ndimension;
			dint id_y = mesh->dim_y + cell*mesh->ndimension;
			dint id_z = mesh->dim_z + cell*mesh->ndimension;

			dfloat X = mesh->C[id_x]*lengthScaleFactor;
			dfloat Y = mesh->C[id_y]*lengthScaleFactor;
			dfloat Z = mesh->C[id_z]*lengthScaleFactor;

			dfloat Phi = phi[cell]*fieldScaleFactor;



			fprintf(fptr,"%9.5e\t%9.5e\t%9.5e\t%9.5e\t",X, Y, Z, Phi);

			fprintf(fptr,"\n");


		}

	}

	//------------------------------------------------------------------!
	//  Volume Vector
	//------------------------------------------------------------------!

	if (fieldType == volumeVector ){

		fprintf(fptr,"VARIABLES = \"X\" \"Y\" \"Z\" \"%sX\" \"%sY\" \"%sZ\"\n",mesh->variableName,mesh->variableName,mesh->variableName);

		fprintf(fptr,"ZONE I=%d, DATAPACKING=POINT \t \n",mesh->ncellsx);

		dint j = J;
		dint i = I;

		for (dint k = 0;k < mesh->ncellsz;k++){

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

			fprintf(fptr,"%9.5e\t%9.5e\t%9.5e\t%9.5e\t%9.5e\t%9.5e\t", X, Y, Z, Phi_X, Phi_Y, Phi_Z);

			fprintf(fptr,"\n");


		}

	}

	fclose(fptr);


}


