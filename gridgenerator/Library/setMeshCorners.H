#include "mesh.h"

void setMeshCorners(mesh_t *mesh);

void setMeshCorners(mesh_t *mesh){

	// BOTTOM BOUNDARY ( Z axis is out of Plane )
	/*        	 Y
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
	 * (X1,Y1,Z1)*---------------------------*(X2,Y2,Z2)-----------> X
	 * 				       SOUTH
	 */


	mesh->x1 = mesh->x_min;   mesh->y1 =  mesh->y_min;   mesh->z1 = mesh->z_min;
	mesh->x2 = mesh->x_max;   mesh->y2 =  mesh->y_min;   mesh->z2 = mesh->z_min;
	mesh->x3 = mesh->x_max;   mesh->y3 =  mesh->y_max;   mesh->z3 = mesh->z_min;
	mesh->x4 = mesh->x_min;   mesh->y4 =  mesh->y_max;   mesh->z4 = mesh->z_min;


	// TOP BOUNDARY ( Z axis is out of Plane )
	/*        	 Y
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
	 * (X5,Y5,Z5)*---------------------------*(X6,Y6,Z6)-----------> X
	 * 				       SOUTH
	 */

	mesh->x5 = mesh->x_min;   mesh->y5 =  mesh->y_min;   mesh->z5 = mesh->z_max;
	mesh->x6 = mesh->x_max;   mesh->y6 =  mesh->y_min;   mesh->z6 = mesh->z_max;
	mesh->x7 = mesh->x_max;   mesh->y7 =  mesh->y_max;   mesh->z7 = mesh->z_max;
	mesh->x8 = mesh->x_min;   mesh->y8 =  mesh->y_max;   mesh->z8 = mesh->z_max;

}
