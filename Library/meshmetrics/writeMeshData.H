#include "mesh.h"

void writeMeshData(dfloat* Phi, int fieldType,mesh_t *mesh);

void writeMeshData(dfloat* Phi, int fieldType,mesh_t *mesh){

	strcpy(mesh->directoryName, "Mesh/");
	mkdir(mesh->directoryName, 0777);
	strcat(mesh->directoryName, mesh->variableName);
	strcpy(mesh->variableName, mesh->directoryName);

	dint ndimension   = getNdimension(fieldType);
	dint ndataPerCell = getNdatapercell(fieldType);

	FILE * fptr;
	strcpy(mesh->fileName, mesh->variableName); strcat(mesh->fileName,  ".txt");

	fptr = fopen (mesh->fileName, "w");

	if(fptr == NULL) {
		printf("Could not create file.\n");
		exit(-1);
	}

	//------------------------------------------------------------------!
	//  Volume Scalar
	//------------------------------------------------------------------!

	if (fieldType == volumeScalar ){

		fprintf(fptr,"VARIABLES = \"X\" \"Y\" \"Z\" \"%s\"\n",mesh->variableName);


		fprintf(fptr,"ZONE I=%d, J=%d, K=%d, DATAPACKING=POINT \t \n",mesh->ncellsx,mesh->ncellsy,mesh->ncellsz);

		for(dint k = 0; k < mesh->ncellsz; k++){
			for (dint j = 0;j < mesh->ncellsy;j++){
				for (dint i = 0;i < mesh->ncellsx;i++){

					dint cell = computeLinearIndex(i,j,k, mesh->ncellsx, mesh->ncellsy,mesh->ncellsz);

					dint id_x = mesh->dim_x + cell*mesh->ndimension;
					dint id_y = mesh->dim_y + cell*mesh->ndimension;
					dint id_z = mesh->dim_z + cell*mesh->ndimension;

					fprintf(fptr,"%9.2e\t%9.2e\t%9.2e\t%9.2e\t",mesh->C[id_x], mesh->C[id_y], mesh->C[id_z],Phi[cell] );

					fprintf(fptr,"\n");


				}
			}
		}
	}

	//------------------------------------------------------------------!
	//  Volume Vector
	//------------------------------------------------------------------!

	if (fieldType == volumeVector ){

		char *variableNameX ;
		char *variableNameY ;
		char *variableNameZ ;

		createCharField(&variableNameX,mesh);
		createCharField(&variableNameY,mesh);
		createCharField(&variableNameZ,mesh);

		strcpy(variableNameX, mesh->variableName);
		strcpy(variableNameY, mesh->variableName);
		strcpy(variableNameZ, mesh->variableName);

		strcat(variableNameX, ".X");
		strcat(variableNameY, ".Y");
		strcat(variableNameZ, ".Z");


		fprintf(fptr,"VARIABLES = \"X\" \"Y\" \"Z\" \"%s\" \"%s\" \"%s\"\n",variableNameX,variableNameY,variableNameZ);

		fprintf(fptr,"ZONE I=%d, J=%d, K=%d, DATAPACKING=POINT \t \n",mesh->ncellsx,mesh->ncellsy,mesh->ncellsz);

		for(dint k = 0; k < mesh->ncellsz; k++){
			for (dint j = 0;j < mesh->ncellsy;j++){
				for (dint i = 0;i < mesh->ncellsx;i++){

					dint cell = computeLinearIndex(i,j,k, mesh->ncellsx, mesh->ncellsy,mesh->ncellsz);

					dint id_x = mesh->dim_x + cell*mesh->ndimension;
					dint id_y = mesh->dim_y + cell*mesh->ndimension;
					dint id_z = mesh->dim_z + cell*mesh->ndimension;

					fprintf(fptr,"%9.2e\t%9.2e\t%9.2e\t%9.2e\t%9.2e\t%9.2e\t",mesh->C[id_x], mesh->C[id_y], mesh->C[id_z],Phi[id_x],Phi[id_y],Phi[id_z] );

					fprintf(fptr,"\n");


				}
			}
		}


		free(variableNameX);
		free(variableNameY);
		free(variableNameZ);
	}

	//------------------------------------------------------------------!
	//  Surface Scalar
	//------------------------------------------------------------------!

	if (fieldType == surfaceScalar ){

		fprintf(fptr,"VARIABLES = \"X\" \"Y\" \"Z\" \"faceAEFB\" \"faceEGHF\" \"faceGCDH\" \"faceCABD\" \"faceFHDB\" \"faceACGE\"\n");

		fprintf(fptr,"ZONE I=%d, J=%d, K=%d, DATAPACKING=POINT \t \n",mesh->ncellsx,mesh->ncellsy,mesh->ncellsz);

		for(dint k = 0; k < mesh->ncellsz; k++){
			for (dint j = 0;j < mesh->ncellsy;j++){
				for (dint i = 0;i < mesh->ncellsx;i++){

					dint cell = computeLinearIndex(i,j,k, mesh->ncellsx, mesh->ncellsy,mesh->ncellsz);

					dint id_x = mesh->dim_x + cell*mesh->ndimension;
					dint id_y = mesh->dim_y + cell*mesh->ndimension;
					dint id_z = mesh->dim_z + cell*mesh->ndimension;

					dint fidAEFB  =  mesh->faceAEFB + cell*mesh->ncellFaces;
					dint fidEGHF  =  mesh->faceEGHF + cell*mesh->ncellFaces;
					dint fidGCDH  =  mesh->faceGCDH + cell*mesh->ncellFaces;
					dint fidCABD  =  mesh->faceCABD + cell*mesh->ncellFaces;
					dint fidFHDB  =  mesh->faceFHDB + cell*mesh->ncellFaces;
					dint fidACGE  =  mesh->faceACGE + cell*mesh->ncellFaces;



					fprintf(fptr,"%9.2e\t%9.2e\t%9.2e\t%9.2e\t%9.2e\t%9.2e\t%9.2e\t%9.2e\t%9.2e\t",mesh->C[id_x], mesh->C[id_y], mesh->C[id_z],
							     Phi[fidAEFB],Phi[fidEGHF],Phi[fidGCDH],Phi[fidCABD],Phi[fidFHDB],Phi[fidACGE] );

					fprintf(fptr,"\n");


				}
			}
		}
	}

	//------------------------------------------------------------------!
	//  Surface Vector
	//------------------------------------------------------------------!

	if (fieldType == surfaceVector ){

		fprintf(fptr,"VARIABLES = \"X\" \"Y\" \"Z\""
				" \"faceAEFB_x\" \"faceAEFB_y\" \"faceAEFB_z\""
				" \"faceEGHF_x\" \"faceEGHF_y\" \"faceEGHF_z\""
				" \"faceGCDH_x\" \"faceGCDH_y\" \"faceGCDH_z\""
				" \"faceCABD_x\" \"faceCABD_y\" \"faceCABD_z\""
				" \"faceFHDB_x\" \"faceFHDB_y\" \"faceFHDB_z\""
				" \"faceACGE_x\" \"faceACGE_y\" \"faceACGE_z\"\n");

		fprintf(fptr,"ZONE I=%d, J=%d, K=%d, DATAPACKING=POINT \t \n",mesh->ncellsx,mesh->ncellsy,mesh->ncellsz);

		for(dint k = 0; k < mesh->ncellsz; k++){
			for (dint j = 0;j < mesh->ncellsy;j++){
				for (dint i = 0;i < mesh->ncellsx;i++){

					dint cell = computeLinearIndex(i,j,k, mesh->ncellsx, mesh->ncellsy,mesh->ncellsz);

					dint id_x = mesh->dim_x + cell*mesh->ndimension;
					dint id_y = mesh->dim_y + cell*mesh->ndimension;
					dint id_z = mesh->dim_z + cell*mesh->ndimension;

					dint fidAEFB  =  mesh->faceAEFB + cell*mesh->ncellFaces;
					dint fidEGHF  =  mesh->faceEGHF + cell*mesh->ncellFaces;
					dint fidGCDH  =  mesh->faceGCDH + cell*mesh->ncellFaces;
					dint fidCABD  =  mesh->faceCABD + cell*mesh->ncellFaces;
					dint fidFHDB  =  mesh->faceFHDB + cell*mesh->ncellFaces;
					dint fidACGE  =  mesh->faceACGE + cell*mesh->ncellFaces;

					dint fidAEFB_x  =  mesh->dim_x + fidAEFB*mesh->ndimension;
					dint fidAEFB_y  =  mesh->dim_y + fidAEFB*mesh->ndimension;
					dint fidAEFB_z  =  mesh->dim_z + fidAEFB*mesh->ndimension;

					dint fidEGHF_x  =  mesh->dim_x + fidEGHF*mesh->ndimension;
					dint fidEGHF_y  =  mesh->dim_y + fidEGHF*mesh->ndimension;
					dint fidEGHF_z  =  mesh->dim_z + fidEGHF*mesh->ndimension;

					dint fidGCDH_x  =  mesh->dim_x + fidGCDH*mesh->ndimension;
					dint fidGCDH_y  =  mesh->dim_y + fidGCDH*mesh->ndimension;
					dint fidGCDH_z  =  mesh->dim_z + fidGCDH*mesh->ndimension;

					dint fidCABD_x  =  mesh->dim_x + fidCABD*mesh->ndimension;
					dint fidCABD_y  =  mesh->dim_y + fidCABD*mesh->ndimension;
					dint fidCABD_z  =  mesh->dim_z + fidCABD*mesh->ndimension;

					dint fidFHDB_x  =  mesh->dim_x + fidFHDB*mesh->ndimension;
					dint fidFHDB_y  =  mesh->dim_y + fidFHDB*mesh->ndimension;
					dint fidFHDB_z  =  mesh->dim_z + fidFHDB*mesh->ndimension;

					dint fidACGE_x  =  mesh->dim_x + fidACGE*mesh->ndimension;
					dint fidACGE_y  =  mesh->dim_y + fidACGE*mesh->ndimension;
					dint fidACGE_z  =  mesh->dim_z + fidACGE*mesh->ndimension;

					fprintf(fptr,"%9.2e\t%9.2e\t%9.2e\t"
							     "%9.2e\t%9.2e\t%9.2e\t"
							     "%9.2e\t%9.2e\t%9.2e\t"
							     "%9.2e\t%9.2e\t%9.2e\t"
							     "%9.2e\t%9.2e\t%9.2e\t"
							     "%9.2e\t%9.2e\t%9.2e\t"
							     "%9.2e\t%9.2e\t%9.2e\t",
								 mesh->C[id_x], mesh->C[id_y], mesh->C[id_z],
							     Phi[fidAEFB_x],Phi[fidAEFB_y],Phi[fidAEFB_z],
								 Phi[fidEGHF_x],Phi[fidEGHF_y],Phi[fidEGHF_z],
								 Phi[fidGCDH_x],Phi[fidGCDH_y],Phi[fidGCDH_z],
								 Phi[fidCABD_x],Phi[fidCABD_y],Phi[fidCABD_z],
								 Phi[fidFHDB_x],Phi[fidFHDB_y],Phi[fidFHDB_z],
								 Phi[fidACGE_x],Phi[fidACGE_y],Phi[fidACGE_z] );

					fprintf(fptr,"\n");


				}
			}
		}
	}

	//------------------------------------------------------------------!
	//  Vertex Scalar
	//------------------------------------------------------------------!

	if (fieldType == vertexScalar ){

		fprintf(fptr,"VARIABLES = \"X\" \"Y\" \"Z\" \"vertA\" \"vertB\" \"vertC\" \"vertD\" \"vertE\" \"vertF\" \"vertG\" \"vertH\"\n");

		fprintf(fptr,"ZONE I=%d, J=%d, K=%d, DATAPACKING=POINT \t \n",mesh->ncellsx,mesh->ncellsy,mesh->ncellsz);

		for(dint k = 0; k < mesh->ncellsz; k++){
			for (dint j = 0;j < mesh->ncellsy;j++){
				for (dint i = 0;i < mesh->ncellsx;i++){

					dint cell = computeLinearIndex(i,j,k, mesh->ncellsx, mesh->ncellsy,mesh->ncellsz);

					dint id_x = mesh->dim_x + cell*mesh->ndimension;
					dint id_y = mesh->dim_y + cell*mesh->ndimension;
					dint id_z = mesh->dim_z + cell*mesh->ndimension;

					dint vidvertA  =  mesh->vertA + cell*mesh->ncellVertices;
					dint vidvertB  =  mesh->vertB + cell*mesh->ncellVertices;
					dint vidvertC  =  mesh->vertC + cell*mesh->ncellVertices;
					dint vidvertD  =  mesh->vertD + cell*mesh->ncellVertices;
					dint vidvertE  =  mesh->vertE + cell*mesh->ncellVertices;
					dint vidvertF  =  mesh->vertF + cell*mesh->ncellVertices;
					dint vidvertG  =  mesh->vertG + cell*mesh->ncellVertices;
					dint vidvertH  =  mesh->vertH + cell*mesh->ncellVertices;


					fprintf(fptr,"%9.2e\t%9.2e\t%9.2e\t"
								 "%9.2e\t%9.2e\t%9.2e\t%9.2e\t"
								 "%9.2e\t%9.2e\t%9.2e\t%9.2e\t",
								 mesh->C[id_x],mesh->C[id_y],mesh->C[id_z],
							     Phi[vidvertA],Phi[vidvertB],Phi[vidvertC],Phi[vidvertD],
								 Phi[vidvertE],Phi[vidvertF],Phi[vidvertG],Phi[vidvertH]);

					fprintf(fptr,"\n");


				}
			}
		}
	}

	//------------------------------------------------------------------!
	//  Vertex Vector
	//------------------------------------------------------------------!

	if (fieldType == vertexVector ){

		fprintf(fptr,"VARIABLES = \"X\" \"Y\" \"Z\""
				" \"vertA_x\" \"vertA_y\" \"vertA_z\""
				" \"vertB_x\" \"vertB_y\" \"vertB_z\""
				" \"vertC_x\" \"vertC_y\" \"vertC_z\""
				" \"vertD_x\" \"vertD_y\" \"vertD_z\""
				" \"vertE_x\" \"vertE_y\" \"vertE_z\""
				" \"vertF_x\" \"vertF_y\" \"vertF_z\""
				" \"vertG_x\" \"vertG_y\" \"vertG_z\""
				" \"vertH_x\" \"vertH_y\" \"vertH_z\"\n");

		fprintf(fptr,"ZONE I=%d, J=%d, K=%d, DATAPACKING=POINT \t \n",mesh->ncellsx,mesh->ncellsy,mesh->ncellsz);

		for(dint k = 0; k < mesh->ncellsz; k++){
			for (dint j = 0;j < mesh->ncellsy;j++){
				for (dint i = 0;i < mesh->ncellsx;i++){

					dint cell = computeLinearIndex(i,j,k, mesh->ncellsx, mesh->ncellsy,mesh->ncellsz);

					dint id_x = mesh->dim_x + cell*mesh->ndimension;
					dint id_y = mesh->dim_y + cell*mesh->ndimension;
					dint id_z = mesh->dim_z + cell*mesh->ndimension;

					dint vidvertA  =  mesh->vertA + cell*mesh->ncellVertices;
					dint vidvertB  =  mesh->vertB + cell*mesh->ncellVertices;
					dint vidvertC  =  mesh->vertC + cell*mesh->ncellVertices;
					dint vidvertD  =  mesh->vertD + cell*mesh->ncellVertices;
					dint vidvertE  =  mesh->vertE + cell*mesh->ncellVertices;
					dint vidvertF  =  mesh->vertF + cell*mesh->ncellVertices;
					dint vidvertG  =  mesh->vertG + cell*mesh->ncellVertices;
					dint vidvertH  =  mesh->vertH + cell*mesh->ncellVertices;

					dint vidvertA_x =  mesh->dim_x + vidvertA*mesh->ndimension;
					dint vidvertA_y =  mesh->dim_y + vidvertA*mesh->ndimension;
					dint vidvertA_z =  mesh->dim_z + vidvertA*mesh->ndimension;

					dint vidvertB_x =  mesh->dim_x + vidvertB*mesh->ndimension;
					dint vidvertB_y =  mesh->dim_y + vidvertB*mesh->ndimension;
					dint vidvertB_z =  mesh->dim_z + vidvertB*mesh->ndimension;

					dint vidvertC_x =  mesh->dim_x + vidvertC*mesh->ndimension;
					dint vidvertC_y =  mesh->dim_y + vidvertC*mesh->ndimension;
					dint vidvertC_z =  mesh->dim_z + vidvertC*mesh->ndimension;

					dint vidvertD_x =  mesh->dim_x + vidvertD*mesh->ndimension;
					dint vidvertD_y =  mesh->dim_y + vidvertD*mesh->ndimension;
					dint vidvertD_z =  mesh->dim_z + vidvertD*mesh->ndimension;

					dint vidvertE_x =  mesh->dim_x + vidvertE*mesh->ndimension;
					dint vidvertE_y =  mesh->dim_y + vidvertE*mesh->ndimension;
					dint vidvertE_z =  mesh->dim_z + vidvertE*mesh->ndimension;

					dint vidvertF_x =  mesh->dim_x + vidvertF*mesh->ndimension;
					dint vidvertF_y =  mesh->dim_y + vidvertF*mesh->ndimension;
					dint vidvertF_z =  mesh->dim_z + vidvertF*mesh->ndimension;

					dint vidvertG_x =  mesh->dim_x + vidvertG*mesh->ndimension;
					dint vidvertG_y =  mesh->dim_y + vidvertG*mesh->ndimension;
					dint vidvertG_z =  mesh->dim_z + vidvertG*mesh->ndimension;

					dint vidvertH_x =  mesh->dim_x + vidvertH*mesh->ndimension;
					dint vidvertH_y =  mesh->dim_y + vidvertH*mesh->ndimension;
					dint vidvertH_z =  mesh->dim_z + vidvertH*mesh->ndimension;


					fprintf(fptr,"%9.2e\t%9.2e\t%9.2e\t"
								 "%9.2e\t%9.2e\t%9.2e\t"
								 "%9.2e\t%9.2e\t%9.2e\t"
								 "%9.2e\t%9.2e\t%9.2e\t"
								 "%9.2e\t%9.2e\t%9.2e\t"
								 "%9.2e\t%9.2e\t%9.2e\t"
								 "%9.2e\t%9.2e\t%9.2e\t"
								 "%9.2e\t%9.2e\t%9.2e\t"
								 "%9.2e\t%9.2e\t%9.2e\t",
								 mesh->C[id_x],mesh->C[id_y],mesh->C[id_z],
							     Phi[vidvertA_x],Phi[vidvertA_y],Phi[vidvertA_z],
								 Phi[vidvertB_x],Phi[vidvertB_y],Phi[vidvertB_z],
								 Phi[vidvertC_x],Phi[vidvertC_y],Phi[vidvertC_z],
								 Phi[vidvertD_x],Phi[vidvertD_y],Phi[vidvertD_z],
								 Phi[vidvertE_x],Phi[vidvertE_y],Phi[vidvertE_z],
								 Phi[vidvertF_x],Phi[vidvertF_y],Phi[vidvertF_z],
								 Phi[vidvertG_x],Phi[vidvertG_y],Phi[vidvertG_z],
								 Phi[vidvertH_x],Phi[vidvertH_y],Phi[vidvertH_z]);

					fprintf(fptr,"\n");


				}
			}
		}
	}




	fclose(fptr);

}
