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


	mesh->x1 = mesh->innerRadius;   mesh->y1 = -mesh->axisSymmetricMeshAngle/2.0;   mesh->z1 = -mesh->totalHeight/2.0;
	mesh->x2 = mesh->outerRadius;   mesh->y2 = -mesh->axisSymmetricMeshAngle/2.0;   mesh->z2 = -mesh->totalHeight/2.0;
	mesh->x3 = mesh->outerRadius;   mesh->y3 =  mesh->axisSymmetricMeshAngle/2.0;   mesh->z3 = -mesh->totalHeight/2.0;
	mesh->x4 = mesh->innerRadius;   mesh->y4 =  mesh->axisSymmetricMeshAngle/2.0;   mesh->z4 = -mesh->totalHeight/2.0;


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

	mesh->x5 = mesh->innerRadius;   mesh->y5 =  -mesh->axisSymmetricMeshAngle/2.0;  mesh->z5 = mesh->totalHeight/2.0;
	mesh->x6 = mesh->outerRadius;   mesh->y6 =  -mesh->axisSymmetricMeshAngle/2.0;  mesh->z6 = mesh->totalHeight/2.0;
	mesh->x7 = mesh->outerRadius;   mesh->y7 =   mesh->axisSymmetricMeshAngle/2.0;  mesh->z7 = mesh->totalHeight/2.0;
	mesh->x8 = mesh->innerRadius;   mesh->y8 =   mesh->axisSymmetricMeshAngle/2.0;  mesh->z8 = mesh->totalHeight/2.0;

}
