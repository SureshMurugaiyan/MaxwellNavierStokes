#include "mesh.h"

//#include "../Library/solvers/mainNavierStokesSolve.H"

#include "../setUpCase/FieldNames.H"

void mainNavierStokesSolve(mesh_t *mesh);

void mainNavierStokesSolve(mesh_t *mesh){

	meshConstructConstantsFields(mesh); //OK

	setConstants(mesh); //OK

	setScaleFactors(mesh); //OK

	meshConstructGravityMetricsFields(mesh); //OK

	setReferenceHeight(mesh); //OK

	meshGenerateCellMetricsGravity(mesh); //OK

	constructSolverFields(mesh); //OK

	solverSetControls(mesh); //OK

	mesh->timeStep 		= mesh->starttimeStep ;
	mesh->time     		= 0;
	mesh->stopCriterion	= false;

	constructSolutionFields(mesh); //OK

	printf("\nStarting time marching loop \n");

	while (mesh->stopCriterion==false){

		if(mesh->timeStep  == 0){

			printf("\n!=====================================================================! \n");

			printf("\nTimeStep: %d  Time: %9.5e seconds  \n", mesh->timeStep,mesh->time);

			printf("\n!=====================================================================! \n");

			setInitialConditions(mesh); //OK

			updateHydroStaticPressure(mesh); //OK

			computeInitialVelocityFlux(mesh); //OK

			computeDensityFaceValues(mesh); //OK

			computeDynamicViscosityFaceValues(mesh); //OK

			computeConductivityFaceValues(mesh); //OK

			computePostProcessingParameters(mesh); //OK

			storeOldTimeValues(mesh); //OK

			writeResults(mesh); //OK

		}

		storeOldTimeValues(mesh); //OK

		computeDeltaTime(mesh); //OK

		mesh->time = mesh->time + mesh->deltatime;

		mesh->timeStep++;

		printf("\n!=====================================================================! \n");

		printf("\nTimeStep: %d  Time: %9.5e seconds  \n", mesh->timeStep,mesh->time);

		printf("\n!=====================================================================! \n");

		printf("\n\t Maximum Courant Number Using Convective Velocity CoU: %9.2e  \n", mesh->CoU);

		printf("\n\t Time Step Size (delta time): %9.2e seconds  \n\n", mesh->deltatime);




		constructAlphaEquation(mesh); //OK

		printf("\n\t Solving alpha Equation \n");

		PetscMatrixAssemblyAndSolveScalarEqn(mesh->alpha,mesh); //OK

		updateBoundaryValues(mesh->alpha, alpha, mesh); //OK

		updateMixtureDensity(mesh); //OK

		updateMixtureViscosity(mesh); //OK

		updateMixtureConductivity(mesh); //OK

		updateHydroStaticPressure(mesh); //OK



		setUpExternalMagneticField(mesh);

		computeExternalMagneticFieldFlux(mesh);

		// Solving Maxwell Equations

		constructPhiEquation(mesh);

		printf("\n\t Solving Phi Equation \n");

		PetscMatrixAssemblyAndSolveScalarEqn(mesh->phi,mesh);

		updateBoundaryValues(mesh->phi, phi, mesh);

		computeResidualPhiEquation(mesh);

		computeCurrentDensityFlux(mesh);

		computeCurrentDensityAtCellCenter(mesh);



		storeReferenceVelocity(mesh); //OK

		computeBodyForces(mesh); //OK

		constructUEquation(mesh); //OK

		relaxUEqn(mesh); //OK

		printf("\n\t Solving U Equation \n");

		PetscMatrixAssemblyAndSolveVectorEqn(mesh->U,mesh); //OK

		updateBoundaryValues(mesh->U, U, mesh); //OK

		computeResidualUEquation(mesh); //OK

		constructrAP(mesh); //OK


		// Solving Navier-Stokes Equations

		mesh->currentNCorrector = 1;

		while ( mesh->currentNCorrector <= mesh->targetNCorrector){

			storeReferencePressure(mesh); //OK

			constructHbyA(mesh); //OK

			constructPEquation(mesh); //OK

			relaxPEqn(mesh); //OK

			printf("\n\t Solving P Equation \n");

			PetscMatrixAssemblyAndSolveScalarEqn(mesh->p,mesh); //OK

			updateBoundaryValues(mesh->p, p, mesh); //OK

			setReferencePressure(mesh);//OK






			computeResidualPEquation(mesh); //OK

			constructGradPbyA(mesh); //OK

			correctU(mesh); //OK

			correctUFluxes(mesh); //OK

			mesh->currentNCorrector++;

		}


		computePostProcessingParameters(mesh); //OK

		writeResults(mesh); //OK

		solverUpdateStoppingCriterion(mesh); //OK


	}





	destructSolutionFields(mesh);

	destructSolverFields(mesh);

	meshDestructGravityMetricsFields(mesh);

	meshDestructConstantsFields(mesh);

}



