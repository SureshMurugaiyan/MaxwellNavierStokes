#include "mesh.h"

void setMeshControls(mesh_t *mesh);

void setMeshControls(mesh_t *mesh){

	mesh->ncellsx = 50;
	mesh->ncellsy = 1;
	mesh->ncellsz = 26;

	mesh->convertCylindricalPolarToCartesian = true;
	mesh->convertSphericalPolarToCartesian   = false;



	mesh->axisSymmetricMeshAngle = 5*PI/180.0;

	mesh->innerRadius = 1.00;
	mesh->outerRadius = 5.00;
	mesh->totalHeight = 2.00;


	// Z axis is out of Plane
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



//
//	printf ( "!--------------------------------- !\n" );
//	printf ( "! setMeshControls                  !\n" );
//	printf ( "!--------------------------------- !\n!" );
//	printf ( "%13s %5d \n!","ncellsx = ",mesh->ncellsx);
//	printf ( "%13s %5d \n!","ncellsy = ",mesh->ncellsy);
//	printf ( "%13s %5d \n!","ncellsz = ",mesh->ncellsz);
//	printf("convertCylindricalPolarToCartesian = %s\n!", mesh->convertCylindricalPolarToCartesian ? "true" : "false");
//	printf("  convertSphericalPolarToCartesian = %s\n!", mesh->convertSphericalPolarToCartesian ? "true" : "false");
//	printf ( "--------------------------------- !\n" );
//
//	printf ( "\n");




}


