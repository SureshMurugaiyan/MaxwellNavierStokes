#include "mesh.h"

//#include "../setUpCase/computePostProcessingParameters.H"

void computePostProcessingParameters(mesh_t *mesh);

void computePostProcessingParameters(mesh_t *mesh){


	//======================================================!
	// Computing Integrated Total Mass                      !
	//======================================================!

	dfloat sumIMass = 0.0;

	for(int cell = 0;cell < mesh->nicells; ++cell){

		dfloat LiquidMass = mesh->alpha[cell]*mesh->rhoL*mesh->V[cell];
		dfloat GasMass =  (1- mesh->alpha[cell])*mesh->rhoG*mesh->V[cell];


		dfloat Mass = LiquidMass + GasMass;

		sumIMass = sumIMass + Mass;

	}

	mesh->iMass = sumIMass;

	printf ( "\n\n\t Integrated Mass  :  %9.10e \n\n",mesh->iMass );

	//======================================================!
	// Computing Integrated  Energy                         !
	//======================================================!

	dfloat sumKE  = 0.0; // Integrating Kinetic Energy

	dfloat sumKEx = 0.0; // Integrating Kinetic Energy - x component
	dfloat sumKEy = 0.0; // Integrating Kinetic Energy - y component
	dfloat sumKEz = 0.0; // Integrating Kinetic Energy - z component

	for(int cell = 0;cell < mesh->nicells; ++cell){

		int idx,idy,idz;

		dfloat x,y,z,u,v,w;

		idx	= mesh->dim_x+cell*mesh->ndimension;
		idy	= mesh->dim_y+cell*mesh->ndimension;
		idz	= mesh->dim_z+cell*mesh->ndimension;

		x = mesh->C[idx];
		y = mesh->C[idy];
		z = mesh->C[idz];

		u = mesh->U[idx];
		v = mesh->U[idy];
		w = mesh->U[idz];

		dfloat KE	= 0.5*(u*u + v*v + w*w)*mesh->rho[cell]*mesh->V[cell];

		dfloat KEx	= 0.5*(u*u)*mesh->rho[cell]*mesh->V[cell];
		dfloat KEy	= 0.5*(v*v)*mesh->rho[cell]*mesh->V[cell];
		dfloat KEz	= 0.5*(w*w)*mesh->rho[cell]*mesh->V[cell];

		sumKE  	 = sumKE  + KE ;
		sumKEx 	 = sumKEx + KEx ;
		sumKEy   = sumKEy + KEy ;
		sumKEz 	 = sumKEz + KEz ;

	}

	mesh->iKE   = sumKE;
	mesh->iKEx  = sumKEx;
	mesh->iKEy  = sumKEy;
	mesh->iKEz  = sumKEz;


	printf ( "\n\n\t Integrated Kinetic Energy  :  %9.10e",mesh->iKE );

	printf ( "\n\n\t Integrated Kinetic Energy x :  %9.10e", mesh->iKEx );
	printf (   "\n\t Integrated Kinetic Energy y :  %9.10e", mesh->iKEy );
	printf (   "\n\t Integrated Kinetic Energy z :  %9.10e \n\n", mesh->iKEz );

	//===============================================================================================================!

	mesh->CPUelapsedtime = clock() - mesh->CPUstarttime;

	mesh->CPUelapsedtimeInSeconds = ((dfloat)mesh->CPUelapsedtime)/CLOCKS_PER_SEC; // in seconds


}
