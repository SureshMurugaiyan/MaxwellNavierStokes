#include "mesh.h"
//#include "../setUpCase/setScaleFactors.H"

void setScaleFactors(mesh_t *mesh);

void setScaleFactors(mesh_t *mesh){

	mesh->L_0		=	1.0;
	mesh->rho_0		=	1.0;
	mesh->g_0		=	1.0;
	mesh->U_0		=	1.0;
	mesh->t_0		=	1.0;
	mesh->p_0		=	1.0;
	mesh->mu_0		=	1.0;


	mesh->sigma_0   =	1.0;
	mesh->Phi_0     =	1.0;
	mesh->B_0 	    =	1.0;
	mesh->J_0 	    =	1.0;

}

