#include "mesh.h"

void generateSouthEastCornerEdge(mesh_t *mesh);

void generateSouthEastCornerEdge(mesh_t *mesh){

	// South East Corner Edge, Line Joining (X2,Y2,Z2) and (X6,Y6,Z6)

	// i = 0                    --> West Boundary
	// i = (mesh->nverticesx-1) --> East Boundary
	// j = 0                    --> South Boundary
	// j = (mesh->nverticesy-1) --> North Boundary
	// k = 0                    --> Bottom Plane
	// k = (mesh->nverticesz-1) --> Top Plane

	// BOTTOM BOUNDARY ( Z axis is out of Plane )
	/*        	 Y(j)
	 *  	  	 ^
	 *  	  	 |
	 * 	         |			NORTH
	 * (X4,Y4,Z4)*--------------------------*(X3,Y3,Z3)
	 *  	     |                          |
	 * 	     W	 |                          |  E
	 *       E   |			                |  A
	 *       S   |                          |  S
	 *       T	 |                          |  T
	 *  	     |                          |
	 * (X1,Y1,Z1)*--------------------------#(X2,Y2,Z2)-----------> X(i)
	 * 				       SOUTH
	 */

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
	 * (X5,Y5,Z5)*--------------------------#(X6,Y6,Z6)-----------> X(i)
	 * 				       SOUTH
	 */

	dint i = (mesh->nverticesx-1);
	dint j = 0;

	for( dint k = 1; k < (mesh->nverticesz-1); k++) { // Skipping the first and last corner points

		dint vid = computeLinearIndex(i,j,k, mesh->nverticesx, mesh->nverticesy,mesh->nverticesz);

		dint vid_x = mesh->dim_x + vid*mesh->ndimension;
		dint vid_y = mesh->dim_y + vid*mesh->ndimension;
		dint vid_z = mesh->dim_z + vid*mesh->ndimension;

		dfloat x = mesh->x_max;            		// Parameterized as a function of Zeta
		dfloat y = mesh->y_min;            		// Parameterized as a function of Zeta
		dfloat z = mesh->z_min + mesh->zeta[k]*(mesh->z_max-mesh->z_min);   // Parameterized as a function of Zeta


		mesh->vertices[vid_x] =  x;
		mesh->vertices[vid_y] =  y;
		mesh->vertices[vid_z] =  z;


	}







//	printf ( "!---------------------------------------------------------------------------------------------------- !\n" );
//	printf ( "! generate generateSouthEastCornerEdge                                                                !\n" );
//	printf ( "!---------------------------------------------------------------------------------------------------- !\n" );
//
//	printVertexVector(mesh->vertices,mesh);

}
