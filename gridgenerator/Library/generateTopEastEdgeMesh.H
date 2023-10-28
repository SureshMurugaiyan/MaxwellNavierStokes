#include "mesh.h"

void generateTopEastEdgeMesh(mesh_t *mesh);

void generateTopEastEdgeMesh(mesh_t *mesh){

	// Top East Edge, Line Joining (X6,Y6,Z6)and (X7,Y7,Z7)

	// i = 0                    --> West Boundary
	// i = (mesh->nverticesx-1) --> East Boundary
	// j = 0                    --> South Boundary
	// j = (mesh->nverticesy-1) --> North Boundary
	// k = 0                    --> Bottom Plane
	// k = (mesh->nverticesz-1) --> Top Plane

	// TOP BOUNDARY ( Z axis is out of Plane )
	/*        	 Y(j)
	 *  	  	 ^
	 *  	  	 |
	 * 	         |			NORTH
	 * (X8,Y8,Z8)*--------------------------*(X7,Y7,Z7)
	 *  	     |                          |
	 * 	     W	 |                          |  E
	 *       E   |			                |  A
	 *       S   |                          |  S
	 *       T	 |                          |  T
	 *  	     |                          |
	 * (X5,Y5,Z5)*--------------------------*(X6,Y6,Z6)-----------> X(i)
	 * 				       SOUTH
	 */

	dint i = mesh->nverticesx-1;
	dint k = mesh->nverticesz-1;

	for (dint j = 1;j < (mesh->nverticesy-1);j++){  // Skipping the first and last corner points

		dint vid = computeLinearIndex(i,j,k, mesh->nverticesx, mesh->nverticesy,mesh->nverticesz);

		dint vid_x = mesh->dim_x + vid*mesh->ndimension;
		dint vid_y = mesh->dim_y + vid*mesh->ndimension;
		dint vid_z = mesh->dim_z + vid*mesh->ndimension;

		dfloat x = mesh->x_max;                // Parameterized as a function of eta
		dfloat y = mesh->y_min + mesh->eta[j]*(mesh->y_max - mesh->y_min);    // Parameterized as a function of eta
		dfloat z = mesh->z_max;                // Parameterized as a function of eta

		mesh->vertices[vid_x] =  x;
		mesh->vertices[vid_y] =  y;
		mesh->vertices[vid_z] =  z;

	}

}
