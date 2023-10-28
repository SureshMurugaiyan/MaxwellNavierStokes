#include "mesh.h"

//#include "../setUpCase/solverSetControls.H"

void solverSetControls(mesh_t *mesh);

void solverSetControls(mesh_t *mesh){

	// Time Loop Controls
	mesh->starttimeStep     	= 0; // start time step
	mesh->timeStepMax			= 25000; //250000
	mesh->finalTime  			= 25;
	mesh->deltatime  			= 1E-3;
	mesh->WriteFrequency		= 100; //100

	mesh->timeStepSizeControl	= initialTimeStep; // choose bt 1.targetCourantNumber 2.initialTimeStep
	mesh->Cot					= 0.5; // target courant number
	mesh->Kl					= 1.05; // under relaxation factor for time step increment
	mesh->lambdat				= 0.10; // limit for overshoots time step

	mesh->divScheme				= Gamma;
	mesh->divSchemeblendFactor	= 0.2;

	mesh->targetNCorrector      = 3;

	mesh->UEqnlambda			= 0.9;
	mesh->PEqnlambda			= 0.8;
	mesh->alphaEqnlambda		= 1.0;


	mesh->includeGravity		= false;
	mesh->includeLorentzForce	= false;


	//=========================================!
	// Notes
	//=========================================!

	/*  divSchemeblendFactor is blending factor for divergence scheme
	 *  divSchemeblendFactor =  0.1, represents 10 % central differencing and 90 % upwind
	 *  divSchemeblendFactor =  0.0, represents pure upwind
	 *  divSchemeblendFactor =  1.0, represents pure central differencing
	 * */

	/*  divSchemeblendFactor is blending factor if "gamma " scheme is selected for divergence scheme
	 *  divSchemeblendFactor should be between 0.1 and 0.5 if Gamma scheme is selected
	 *  divSchemeblendFactor = 0.5 gives more numerical diffusion
	 *  divSchemeblendFactor = 0.1 gives sharp profile
	 *  Keep divSchemeblendFactor close to 0.5, if convergence problem occurs
	 * */

	//=========================================!
	// printing
	//=========================================!

	//	printf ( "!--------------------------------------------- !\n" );
	//	printf ( "! Solver Controls                              !\n" );
	//	printf ( "!--------------------------------------------- !\n!" );
	//
	//	printf ( " start timeStep   : %9d   \n!", mesh->starttimeStep);
	//
	//	printf ( " timeStepMax      : %9d   \n!", mesh->timeStepMax);
	//
	//	printf ( " finalTime        : %9.2e \n!", mesh->finalTime);
	//
	//	printf ( " deltatime        : %9.2e \n!", mesh->deltatime);
	//
	//	printf ( " WriteFrequency   : %9d   \n!", mesh->WriteFrequency);
	//	printf ( "--------------------------------------------- !\n" );

}


