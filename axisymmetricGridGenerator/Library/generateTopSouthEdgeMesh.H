#include "mesh.h"

void generateTopSouthEdgeMesh(mesh_t *mesh);

void generateTopSouthEdgeMesh(mesh_t *mesh){


	// Top South Edge, Line Joining (X5,Y5,Z5) and (X6,Y6,Z6)

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


	dint j = 0;
	dint k = mesh->nverticesz-1;

	for (dint i = 1;i < (mesh->nverticesx-1);i++){ // Skipping the first and last corner points

		dint vid = computeLinearIndex(i,j,k, mesh->nverticesx, mesh->nverticesy,mesh->nverticesz);

		dint vid_x = mesh->dim_x + vid*mesh->ndimension;
		dint vid_y = mesh->dim_y + vid*mesh->ndimension;
		dint vid_z = mesh->dim_z + vid*mesh->ndimension;

		dfloat x =  mesh->innerRadius + mesh->xi[i]*(mesh->outerRadius - mesh->innerRadius); // Parameterized as a function of xi
		dfloat y = -mesh->axisSymmetricMeshAngle/2.0; 										 // Parameterized as a function of xi
		dfloat z =  mesh->totalHeight/2.0;            										 // Parameterized as a function of xi

		mesh->vertices[vid_x] =  x;
		mesh->vertices[vid_y] =  y;
		mesh->vertices[vid_z] =  z;


	}


}
