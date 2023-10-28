#include "mesh.h"

void convertSphericalPolarToCartesian(mesh_t *mesh);

void convertSphericalPolarToCartesian(mesh_t *mesh){

	if(mesh->convertSphericalPolarToCartesian == false) {

		return;
	}


	for(int vertices=0;vertices<mesh->nvertices;vertices++){

		dint dim_r = 0; // r
		dint dim_t = 1; // theta
		dint dim_p = 2; // phi

		dint r_id  = dim_r + vertices*mesh->ndimension;
		dint t_id  = dim_t + vertices*mesh->ndimension;
		dint p_id  = dim_p + vertices*mesh->ndimension;

		dfloat r   = mesh->vertices[r_id];
		dfloat t   = mesh->vertices[t_id];
		dfloat p   = mesh->vertices[p_id];

		mesh->vertices[r_id] = r*cos(t)*sin(p);
		mesh->vertices[t_id] = r*sin(t)*sin(p);
		mesh->vertices[p_id] = r*cos(p);

	}

//
//	printf ( "!---------------------------------------------------------------------------------------------------- !\n" );
//	printf ( "! Spherical to Cartesian Conversion                                                                   !\n" );
//	printf ( "!---------------------------------------------------------------------------------------------------- !\n" );
//
//	printVertexVector(mesh->vertices,mesh);
}
