#include "mesh.h"

void meshSettings(mesh_t *mesh);

void meshSettings(mesh_t *mesh){

	mesh->SMALL 	= 1.0e-15;
	mesh->VSMALL 	= 1.0e-150;
	mesh->VVSMALL 	= 1.0e-300;

	mesh->millisecond 	= 1E-3;
	mesh->microsecond 	= 1E-6;
	mesh->nanosecond 	= 1E-9;
	mesh->picosecond 	= 1E-12;
	mesh->femtosecond 	= 1E-15;

	mesh->minute = 60;
	mesh->hour 	 = 60*mesh->minute;
	mesh->day 	 = 24*mesh->hour;

	mesh->PI = 3.1415926535;
	mesh->permeabilityofFreeSpace = 4*mesh->PI*(1E-7);


	/*                  ^ Z(k)
	 *                  .
	 *                  .
     *                  .
     *                 B*--------------------------*D
	 *				   *.                         *|
     *                * .                        * |
	 *      	     *  .                       *  |
	 *  	  	    *   .                      *   |
	 *  	  	   *    .                     *    |
	 * 	          *		.	                 *     |
	 *          F*--------------------------*H     |
	 *  	     |      .                   |      |
	 * 	     	 |      .                   |      |
	 *           |	   A. . . . . . . . . . | . . .*C-----------> Y(j)
	 *           |     .                    |     *
	 *       	 |   .                      |   *
	 *  	     | .                        | *
	 *          E*--------------------------*G
	 * 		  .
	        .
	      .
	     X(i)
	*/


	/*                  ^ Z(k)
	 *                  .
	 *                  .
     *                  .
     *                1 *--------------------------*
	 *				   *.                         *|3
     *                * .                        * |
	 *      	     *  .                       *  |
	 *  	  	    *   .                      *   |
	 *  	  	   *    .                     *    |
	 * 	          *		.	                 *     |
	 *         5 *--------------------------* 7    |
	 *  	     |      .                   |      |
	 * 	     	 |      .                   |      |2
	 *           |	   0. . . . . . . . . . | .    *-----------> Y(j)
	 *           |     .                    |     *
	 *       	 |   .                      |   *
	 *  	     | .                        | *
	 *         4 *--------------------------*6
	 * 		  .
	        .
	      .
	     X(i)
	*/


	mesh->vertA = 0;
	mesh->vertB = 1;
	mesh->vertC = 2;
	mesh->vertD = 3;
	mesh->vertE = 4;
	mesh->vertF = 5;
	mesh->vertG = 6;
	mesh->vertH = 7;

	mesh->faceAEFB = 0;
	mesh->faceEGHF = 1;
	mesh->faceGCDH = 2;
	mesh->faceCABD = 3;
	mesh->faceFHDB = 4;
	mesh->faceACGE = 5;

	mesh->leftFace 	 = 0; // faceAEFB
	mesh->frontFace  = 1; // faceEGHF
	mesh->rightFace  = 2; // faceGCDH
	mesh->backFace   = 3; // faceCABD
	mesh->topFace    = 4; // faceFHDB
	mesh->bottomFace = 5; // faceACGE



//	printf ( "!--------------------------------- !\n" );
//	printf ( "! Mesh Settings                    !\n" );
//	printf ( "!--------------------------------- !\n!" );
//	printf ( "%13s %9.2e \n!","SMALL   = ",mesh->SMALL);
//	printf ( "%13s %9.2e \n!","VSMALL  = ",mesh->VSMALL);
//	printf ( "%13s %9.2e \n!","VVSMALL = ",mesh->VVSMALL);
//	printf ( "--------------------------------- !\n!" );
//
//	printf ( "%13s %9.2e \n!","millisecond = ",mesh->millisecond);
//	printf ( "%13s %9.2e \n!","microsecond = ",mesh->microsecond);
//	printf ( "%13s %9.2e \n!","nanosecond  = ",mesh->nanosecond);
//	printf ( "%13s %9.2e \n!","picosecond  = ",mesh->picosecond);
//	printf ( "%13s %9.2e \n!","femtosecond = ",mesh->femtosecond);
//	printf ( "--------------------------------- !\n!" );
//
//	printf ( "%13s %9.2f \n!","minute = ",mesh->minute);
//	printf ( "%13s %9.2f \n!","hour   = ",mesh->hour);
//	printf ( "%13s %9.2f \n!","day    = ",mesh->day);
//	printf ( "--------------------------------- !\n!" );
//
//	printf ( "%30s %9.10f \n!","PI = ",mesh->PI);
//	printf ( "%30s %9.2e \n!","permeabilityofFreeSpace = ",mesh->permeabilityofFreeSpace);
//	printf ( "--------------------------------- !\n" );

}
