#include "mesh.h"

void setMeshControls(mesh_t *mesh);

void setMeshControls(mesh_t *mesh){

	mesh->ncellsx = 49;
	mesh->ncellsy = 1;
	mesh->ncellsz = 49;

	mesh->x_min =  0.0;
	mesh->x_max =  1.0;

	mesh->y_min = -0.025;
	mesh->y_max =  0.025;

	mesh->z_min =  0.0;
	mesh->z_max =  1.0;


	mesh->convertCylindricalPolarToCartesian = false;
	mesh->convertSphericalPolarToCartesian   = false;


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


