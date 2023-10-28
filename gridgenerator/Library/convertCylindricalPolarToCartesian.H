#include "mesh.h"

void convertCylindricalPolarToCartesian(mesh_t *mesh);

void convertCylindricalPolarToCartesian(mesh_t *mesh){

	if(mesh->convertCylindricalPolarToCartesian == false) {

		return;
	}


	for(int vertices=0;vertices<mesh->nvertices;vertices++){

		dint dim_r = 0; // r
		dint dim_t = 1; // theta
		dint dim_z = 2; // z

		dint r_id  = dim_r + vertices*mesh->ndimension;
		dint t_id  = dim_t + vertices*mesh->ndimension;
		dint z_id  = dim_z + vertices*mesh->ndimension;

		dfloat r   = mesh->vertices[r_id];
		dfloat t   = mesh->vertices[t_id];
		dfloat z   = mesh->vertices[z_id];

		mesh->vertices[r_id] = r*cos(t);
		mesh->vertices[t_id] = r*sin(t);
		mesh->vertices[z_id] = z; // Unchanged, can comment this line. Kept just for illustration

	}

//
//	printf ( "!---------------------------------------------------------------------------------------------------- !\n" );
//	printf ( "! cylindrical to Cartesian Conversion                                                                 !\n" );
//	printf ( "!---------------------------------------------------------------------------------------------------- !\n" );
//
//	printVertexVector(mesh->vertices,mesh);
}
