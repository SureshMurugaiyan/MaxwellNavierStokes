#include "mesh.h"

void calculateIndices(mesh_t *mesh);

void calculateIndices(mesh_t *mesh){

	// BOTTOM BOUNDARY ( Z axis is out of Plane )
	/*         	  Y(j)
	 *  	  	  ^
	 *  	  	  |
	 * 	          |			NORTH
	 *(0,JMAX-1,0)*--------------------------*(IMAX-1,JMAX-1,0)
	 *  	      |                          |
	 * 	      W   |                          |  E
	 *        E   |			                 |  A
	 *        S   |                          |  S
	 *        T	  |                          |  T
	 *  	      |                          |
	 *    (0,0,0)*---------------------------*(IMAX-1,0,0)-----------> X (i)
	 * 				       SOUTH
	 */

	mesh->dim_x = 0;
	mesh->dim_y = 1;
	mesh->dim_z = 2;

	dfloat i,j,k;

	//=============================================================================================================!
	// North West Corner - bottom
	i = 0;
	j = (mesh->nverticesy-1);
	k = 0;
	mesh->northwestbottomVertexID = computeLinearIndex(i,j,k, mesh->nverticesx, mesh->nverticesy,mesh->nverticesz);

	// South West Corner - bottom
	i = 0;
	j = 0;
	k = 0;
	mesh->southwestbottomVertexID = computeLinearIndex(i,j,k, mesh->nverticesx, mesh->nverticesy,mesh->nverticesz);

	// South East Corner - bottom
	i = (mesh->nverticesx-1);
	j = 0;
	k = 0;
	mesh->southeastbottomVertexID = computeLinearIndex(i,j,k, mesh->nverticesx, mesh->nverticesy,mesh->nverticesz);

	// North East Corner - bottom
	i = (mesh->nverticesx-1);
	j = (mesh->nverticesy-1);
	k = 0;
	mesh->northeastbottomVertexID = computeLinearIndex(i,j,k, mesh->nverticesx, mesh->nverticesy,mesh->nverticesz);


	// TOP BOUNDARY ( Z axis is out of Plane )
	/*         	     Y(j)
	 *  	  	      ^
	 *  	  	      |
	 * 	              |			NORTH
	 *(0,JMAX-1,KMAX-1)*--------------------------*(IMAX-1,JMAX-1,KMAX-1)
	 *  	          |                           |
	 * 	          W   |                           |  E
	 *            E   |			                  |  A
	 *            S   |                           |  S
	 *            T	  |                           |  T
	 *  	          |                           |
	 *    (0,0,KMAX-1)*---------------------------*(IMAX-1,0,KMAX-1)-----------> X (i)
	 * 				       SOUTH
	 */



	//=============================================================================================================!
	// North West Corner - top
	i = 0;
	j = (mesh->nverticesy-1);
	k = (mesh->nverticesz-1);
	mesh->northwesttopVertexID = computeLinearIndex(i,j,k, mesh->nverticesx, mesh->nverticesy,mesh->nverticesz);

	// South West Corner - top
	i = 0;
	j = 0;
	k = (mesh->nverticesz-1);
	mesh->southwesttopVertexID = computeLinearIndex(i,j,k, mesh->nverticesx, mesh->nverticesy,mesh->nverticesz);

	// South East Corner - top
	i = (mesh->nverticesx-1);
	j = 0;
	k = (mesh->nverticesz-1);
	mesh->southeasttopVertexID = computeLinearIndex(i,j,k, mesh->nverticesx, mesh->nverticesy,mesh->nverticesz);

	// North East Corner - top
	i = (mesh->nverticesx-1);
	j = (mesh->nverticesy-1);
	k = (mesh->nverticesz-1);
	mesh->northeasttopVertexID = computeLinearIndex(i,j,k, mesh->nverticesx, mesh->nverticesy,mesh->nverticesz);

	//==============================================================================================!
	mesh->northwestbottomVertexID_x = mesh->dim_x + mesh->northwestbottomVertexID*mesh->ndimension ;
	mesh->northwestbottomVertexID_y = mesh->dim_y + mesh->northwestbottomVertexID*mesh->ndimension ;
	mesh->northwestbottomVertexID_z = mesh->dim_z + mesh->northwestbottomVertexID*mesh->ndimension ;

	mesh->southwestbottomVertexID_x = mesh->dim_x + mesh->southwestbottomVertexID*mesh->ndimension ;
	mesh->southwestbottomVertexID_y = mesh->dim_y + mesh->southwestbottomVertexID*mesh->ndimension ;
	mesh->southwestbottomVertexID_z = mesh->dim_z + mesh->southwestbottomVertexID*mesh->ndimension ;

	mesh->southeastbottomVertexID_x = mesh->dim_x + mesh->southeastbottomVertexID*mesh->ndimension ;
	mesh->southeastbottomVertexID_y = mesh->dim_y + mesh->southeastbottomVertexID*mesh->ndimension ;
	mesh->southeastbottomVertexID_z = mesh->dim_z + mesh->southeastbottomVertexID*mesh->ndimension ;

	mesh->northeastbottomVertexID_x = mesh->dim_x + mesh->northeastbottomVertexID*mesh->ndimension ;
	mesh->northeastbottomVertexID_y = mesh->dim_y + mesh->northeastbottomVertexID*mesh->ndimension ;
	mesh->northeastbottomVertexID_z = mesh->dim_z + mesh->northeastbottomVertexID*mesh->ndimension ;

	//==============================================================================================!
	mesh->northwesttopVertexID_x = mesh->dim_x + mesh->northwesttopVertexID*mesh->ndimension ;
	mesh->northwesttopVertexID_y = mesh->dim_y + mesh->northwesttopVertexID*mesh->ndimension ;
	mesh->northwesttopVertexID_z = mesh->dim_z + mesh->northwesttopVertexID*mesh->ndimension ;

	mesh->southwesttopVertexID_x = mesh->dim_x + mesh->southwesttopVertexID*mesh->ndimension ;
	mesh->southwesttopVertexID_y = mesh->dim_y + mesh->southwesttopVertexID*mesh->ndimension ;
	mesh->southwesttopVertexID_z = mesh->dim_z + mesh->southwesttopVertexID*mesh->ndimension ;

	mesh->southeasttopVertexID_x = mesh->dim_x + mesh->southeasttopVertexID*mesh->ndimension ;
	mesh->southeasttopVertexID_y = mesh->dim_y + mesh->southeasttopVertexID*mesh->ndimension ;
	mesh->southeasttopVertexID_z = mesh->dim_z + mesh->southeasttopVertexID*mesh->ndimension ;

	mesh->northeasttopVertexID_x = mesh->dim_x + mesh->northeasttopVertexID*mesh->ndimension ;
	mesh->northeasttopVertexID_y = mesh->dim_y + mesh->northeasttopVertexID*mesh->ndimension ;
	mesh->northeasttopVertexID_z = mesh->dim_z + mesh->northeasttopVertexID*mesh->ndimension ;
	//==============================================================================================!


//
//	printf ( "!--------------------------------- !\n" );
//	printf ( "! Mesh Corner Indices - bottom      !\n" );
//	printf ( "!--------------------------------- !\n!" );
//	printf ( "%13s %5d \n!","northwestVertexID = ",mesh->northwestbottomVertexID);
//	printf ( "%13s %5d \n!","southwestVertexID = ",mesh->southwestbottomVertexID);
//	printf ( "%13s %5d \n!","southeastVertexID = ",mesh->southeastbottomVertexID);
//	printf ( "%13s %5d \n!","northeastVertexID = ",mesh->northeastbottomVertexID);
//	printf ( "--------------------------------- !\n" );
//
//	printf ( "!--------------------------------- !\n" );
//	printf ( "! Mesh Corner Indices - bottom     !\n" );
//	printf ( "!--------------------------------- !\n!" );
//	printf ( "%13s %5d \n!","northwestVertexID,X = ",mesh->northwestbottomVertexID_x);
//	printf ( "%13s %5d \n!","northwestVertexID,Y = ",mesh->northwestbottomVertexID_y);
//	printf ( "%13s %5d \n!","northwestVertexID,Z = ",mesh->northwestbottomVertexID_z);
//	printf ( "--------------------------------- !\n!" );
//	printf ( "%13s %5d \n!","southwestVertexID,X = ",mesh->southwestbottomVertexID_x);
//	printf ( "%13s %5d \n!","southwestVertexID,Y = ",mesh->southwestbottomVertexID_y);
//	printf ( "%13s %5d \n!","southwestVertexID,Z = ",mesh->southwestbottomVertexID_z);
//	printf ( "--------------------------------- !\n!" );
//	printf ( "%13s %5d \n!","southeastVertexID,X = ",mesh->southeastbottomVertexID_x);
//	printf ( "%13s %5d \n!","southeastVertexID,Y = ",mesh->southeastbottomVertexID_y);
//	printf ( "%13s %5d \n!","southeastVertexID,Z = ",mesh->southeastbottomVertexID_z);
//	printf ( "--------------------------------- !\n!" );
//	printf ( "%13s %5d \n!","northeastVertexID,X = ",mesh->northeastbottomVertexID_x);
//	printf ( "%13s %5d \n!","northeastVertexID,Y = ",mesh->northeastbottomVertexID_y);
//	printf ( "%13s %5d \n!","northeastVertexID,Z = ",mesh->northeastbottomVertexID_z);
//	printf ( "--------------------------------- !\n" );
//
//
//	printf ( "!--------------------------------- !\n" );
//	printf ( "! Mesh Corner Indices - top        !\n" );
//	printf ( "!--------------------------------- !\n!" );
//	printf ( "%13s %5d \n!","northwestVertexID = ",mesh->northwesttopVertexID);
//	printf ( "%13s %5d \n!","southwestVertexID = ",mesh->southwesttopVertexID);
//	printf ( "%13s %5d \n!","southeastVertexID = ",mesh->southeasttopVertexID);
//	printf ( "%13s %5d \n!","northeastVertexID = ",mesh->northeasttopVertexID);
//	printf ( "--------------------------------- !\n" );
//
//	printf ( "!--------------------------------- !\n" );
//	printf ( "! Mesh Corner Indices - top        !\n" );
//	printf ( "!--------------------------------- !\n!" );
//	printf ( "%13s %5d \n!","northwestVertexID,X = ",mesh->northwesttopVertexID_x);
//	printf ( "%13s %5d \n!","northwestVertexID,Y = ",mesh->northwesttopVertexID_y);
//	printf ( "%13s %5d \n!","northwestVertexID,Z = ",mesh->northwesttopVertexID_z);
//	printf ( "--------------------------------- !\n!" );
//	printf ( "%13s %5d \n!","southwestVertexID,X = ",mesh->southwesttopVertexID_x);
//	printf ( "%13s %5d \n!","southwestVertexID,Y = ",mesh->southwesttopVertexID_y);
//	printf ( "%13s %5d \n!","southwestVertexID,Z = ",mesh->southwesttopVertexID_z);
//	printf ( "--------------------------------- !\n!" );
//	printf ( "%13s %5d \n!","southeastVertexID,X = ",mesh->southeasttopVertexID_x);
//	printf ( "%13s %5d \n!","southeastVertexID,Y = ",mesh->southeasttopVertexID_y);
//	printf ( "%13s %5d \n!","southeastVertexID,Z = ",mesh->southeasttopVertexID_z);
//	printf ( "--------------------------------- !\n!" );
//	printf ( "%13s %5d \n!","northeastVertexID,X = ",mesh->northeasttopVertexID_x);
//	printf ( "%13s %5d \n!","northeastVertexID,Y = ",mesh->northeasttopVertexID_y);
//	printf ( "%13s %5d \n!","northeastVertexID,Z = ",mesh->northeasttopVertexID_z);
//	printf ( "--------------------------------- !\n" );

}


