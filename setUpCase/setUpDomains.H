#include "mesh.h"
#include "DomainNames.H"

void setUpDomains(mesh_t *mesh);

void setUpDomainsUnitTest(mesh_t *mesh);

void setUpDomains(mesh_t *mesh){

	// Fill the cell center values
	for(int cell = 0;cell < mesh->ncells; ++cell) {

		dint idx	= mesh->dim_x + cell*mesh->ndimension;
		dint idy	= mesh->dim_y + cell*mesh->ndimension;
		dint idz	= mesh->dim_z + cell*mesh->ndimension;

		dfloat x = mesh->C[idx];
		dfloat y = mesh->C[idy];
		dfloat z = mesh->C[idz];


		// Lid Driven Cavity

 		if( z > 0.5){

 			mesh->domain[cell]	= Gas;

 		} else {

 			mesh->domain[cell]	= Liquid;

 		}

//		// Two Phase Lid Driven Cavity
//
// 		if( z < 0.0){
//
// 			mesh->domain[cell]	= Gas;
//
// 		} else {
//
// 			mesh->domain[cell]	= Liquid;
//
// 		}
//
//		// Dam Break Problem
//
//		dfloat xdam =  0.0;
//		dfloat ydam =  0.0;
//
//
// 		if(( x < xdam) && (z < ydam)){
//
// 			mesh->domain[cell]	= Liquid;
//
// 		} else {
//
// 			mesh->domain[cell]	= Gas;
//
// 		}

		// Two Phase MHD

//		if( z > 0){
//
//			mesh->domain[cell]	= Gas;
//
//		} else {
//
//			mesh->domain[cell]	= Liquid;
//
//		}


	}


//	strcpy(mesh->variableName, "domain");
//	printVolumeScalarInt(mesh->domain, mesh);


}
