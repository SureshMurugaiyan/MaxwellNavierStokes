#include "mesh.h"

void assignMeshCorners(mesh_t *mesh);

void assignMeshCorners(mesh_t *mesh){


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

	// Bottom Corners
	mesh->vertices[mesh->northwestbottomVertexID_x] = mesh->x4;// North west Corner
	mesh->vertices[mesh->northwestbottomVertexID_y] = mesh->y4;
	mesh->vertices[mesh->northwestbottomVertexID_z] = mesh->z4;

	mesh->vertices[mesh->southwestbottomVertexID_x] = mesh->x1;// South west Corner
	mesh->vertices[mesh->southwestbottomVertexID_y] = mesh->y1;
	mesh->vertices[mesh->southwestbottomVertexID_z] = mesh->z1;

	mesh->vertices[mesh->southeastbottomVertexID_x] = mesh->x2;// South east Corner
	mesh->vertices[mesh->southeastbottomVertexID_y] = mesh->y2;
	mesh->vertices[mesh->southeastbottomVertexID_z] = mesh->z2;

	mesh->vertices[mesh->northeastbottomVertexID_x] = mesh->x3;// North east Corner
	mesh->vertices[mesh->northeastbottomVertexID_y] = mesh->y3;
	mesh->vertices[mesh->northeastbottomVertexID_z] = mesh->z3;


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

	// Top Corners
	mesh->vertices[mesh->northwesttopVertexID_x] = mesh->x8;// North west Corner
	mesh->vertices[mesh->northwesttopVertexID_y] = mesh->y8;
	mesh->vertices[mesh->northwesttopVertexID_z] = mesh->z8;

	mesh->vertices[mesh->southwesttopVertexID_x] = mesh->x5;// South west Corner
	mesh->vertices[mesh->southwesttopVertexID_y] = mesh->y5;
	mesh->vertices[mesh->southwesttopVertexID_z] = mesh->z5;

	mesh->vertices[mesh->southeasttopVertexID_x] = mesh->x6;// South east Corner
	mesh->vertices[mesh->southeasttopVertexID_y] = mesh->y6;
	mesh->vertices[mesh->southeasttopVertexID_z] = mesh->z6;

	mesh->vertices[mesh->northeasttopVertexID_x] = mesh->x7;// North east Corner
	mesh->vertices[mesh->northeasttopVertexID_y] = mesh->y7;
	mesh->vertices[mesh->northeasttopVertexID_z] = mesh->z7;



//	printf ( "!--------------------------------------------- !\n" );
//	printf ( "! Assigning Mesh corners to global Mesh        !\n" );
//	printf ( "!--------------------------------------------- !\n!" );
//
//	printVertexVector(mesh->vertices,mesh);

}
