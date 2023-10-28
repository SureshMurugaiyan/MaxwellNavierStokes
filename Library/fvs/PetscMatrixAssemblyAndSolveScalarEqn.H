#include "mesh.h"
#include <petscmat.h>
#include "petscksp.h"

//int PetscMatrixAssemblyAndSolveScalarEqn(dfloat* alpha,mesh_t *mesh);

int PetscMatrixAssemblyAndSolveScalarEqn(dfloat* alpha,mesh_t *mesh);

int PetscMatrixAssemblyAndSolveScalarEqn(dfloat* alpha,mesh_t *mesh){

	Mat A; // Linear System Matrix
	Vec X; // Solution Vector
	Vec B; // RHS Source vector
	PetscViewer viewer;  // Viewing the data
	KSP ksp; // Linear solver context
	PC pc;  // Preconditioner context
	PetscErrorCode ierr;
	PetscInt its;

	PetscInt m = mesh->nicells; // number of rows in matrix A
	PetscInt n = mesh->nicells; // number of columns in matrix A
	PetscScalar zero = 0.0;



	// Create Matrix
	ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
	ierr = MatSetType(A,MATSEQAIJ);CHKERRQ(ierr);
	ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m,n);CHKERRQ(ierr);
	ierr = MatSetFromOptions(A);CHKERRQ(ierr);
	ierr = MatSetUp(A);CHKERRQ(ierr);
	ierr = MatZeroEntries(A);CHKERRQ(ierr);


	// Create Solution Vector
	ierr = VecCreate(PETSC_COMM_WORLD,&X);CHKERRQ(ierr);
	ierr = PetscObjectSetName((PetscObject) X, "Solution");CHKERRQ(ierr);
	ierr = VecSetSizes(X,PETSC_DECIDE,n);CHKERRQ(ierr);
	ierr = VecSetFromOptions(X);CHKERRQ(ierr);
	ierr = VecZeroEntries(X);CHKERRQ(ierr);

	// Create RHS Vector
	ierr = VecCreate(PETSC_COMM_WORLD,&B);CHKERRQ(ierr);
	ierr = PetscObjectSetName((PetscObject) B, "RHS Vector");CHKERRQ(ierr);
	ierr = VecSetSizes(B,PETSC_DECIDE,n);CHKERRQ(ierr);
	ierr = VecSetFromOptions(B);CHKERRQ(ierr);
	ierr = VecZeroEntries(B);CHKERRQ(ierr);

	// Matrix Assemble Values for 3D Grid

	for( dint i = 0; i < mesh->ncellsx; i++){
		for( dint j = 0; j < mesh->ncellsy; j++){
			for( dint k = 0; k < mesh->ncellsz; k++){

				dint cell = computeLinearIndex(i,j,k, mesh->ncellsx,mesh->ncellsy,mesh->ncellsz);

				PetscInt row_id = cell;
				PetscInt col_id = cell;
				PetscScalar Avalue = mesh->aP[cell];

				ierr = MatSetValues(A,1,&row_id,1,&col_id,&Avalue,INSERT_VALUES);CHKERRQ(ierr);

			}
		}
	}

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


	// faceAEFB

	for( dint i = 0; i < mesh->ncellsx; i++){
		for( dint j = 1; j < mesh->ncellsy; j++){
			for( dint k = 0; k < mesh->ncellsz; k++){

				dint cell = computeLinearIndex(i,j,k, mesh->ncellsx,mesh->ncellsy,mesh->ncellsz);
				int fid   = mesh->faceAEFB + cell*mesh->ncellFaces;

				PetscInt row_id    = cell;
				PetscInt col_id    = mesh->faceToNeighborCell[fid];
				PetscScalar Avalue = mesh->aN[fid];

				ierr = MatSetValues(A,1,&row_id,1,&col_id,&Avalue,INSERT_VALUES);CHKERRQ(ierr);

			}
		}
	}


	// faceEGHF

	for( dint i = 0; i < (mesh->ncellsx-1); i++){
		for( dint j = 0; j < mesh->ncellsy; j++){
			for( dint k = 0; k < mesh->ncellsz; k++){

				dint cell = computeLinearIndex(i,j,k, mesh->ncellsx,mesh->ncellsy,mesh->ncellsz);
				int fid   = mesh->faceEGHF + cell*mesh->ncellFaces;

				PetscInt row_id    = cell;
				PetscInt col_id    = mesh->faceToNeighborCell[fid];
				PetscScalar Avalue = mesh->aN[fid];

				ierr = MatSetValues(A,1,&row_id,1,&col_id,&Avalue,INSERT_VALUES);CHKERRQ(ierr);

			}
		}
	}


	// faceGCDH

	for( dint i = 0; i < mesh->ncellsx; i++){
		for( dint j = 0; j < (mesh->ncellsy-1); j++){
			for( dint k = 0; k < mesh->ncellsz; k++){

				dint cell = computeLinearIndex(i,j,k, mesh->ncellsx,mesh->ncellsy,mesh->ncellsz);
				int fid   = mesh->faceGCDH + cell*mesh->ncellFaces;

				PetscInt row_id    = cell;
				PetscInt col_id    = mesh->faceToNeighborCell[fid];
				PetscScalar Avalue = mesh->aN[fid];

				ierr = MatSetValues(A,1,&row_id,1,&col_id,&Avalue,INSERT_VALUES);CHKERRQ(ierr);

			}
		}
	}


	// faceCABD

	for( dint i = 1; i < mesh->ncellsx; i++){
		for( dint j = 0; j < mesh->ncellsy; j++){
			for( dint k = 0; k < mesh->ncellsz; k++){

				dint cell = computeLinearIndex(i,j,k, mesh->ncellsx,mesh->ncellsy,mesh->ncellsz);
				int fid   = mesh->faceCABD + cell*mesh->ncellFaces;

				PetscInt row_id    = cell;
				PetscInt col_id    = mesh->faceToNeighborCell[fid];
				PetscScalar Avalue = mesh->aN[fid];

				ierr = MatSetValues(A,1,&row_id,1,&col_id,&Avalue,INSERT_VALUES);CHKERRQ(ierr);

			}
		}
	}


	// faceFHDB

	for( dint i = 0; i < mesh->ncellsx; i++){
		for( dint j = 0; j < mesh->ncellsy; j++){
			for( dint k = 0; k < (mesh->ncellsz-1); k++){

				dint cell = computeLinearIndex(i,j,k, mesh->ncellsx,mesh->ncellsy,mesh->ncellsz);
				int fid   = mesh->faceFHDB + cell*mesh->ncellFaces;

				PetscInt row_id    = cell;
				PetscInt col_id    = mesh->faceToNeighborCell[fid];
				PetscScalar Avalue = mesh->aN[fid];

				ierr = MatSetValues(A,1,&row_id,1,&col_id,&Avalue,INSERT_VALUES);CHKERRQ(ierr);

			}
		}
	}


	// faceACGE

	for( dint i = 0; i < mesh->ncellsx; i++){
		for( dint j = 0; j < mesh->ncellsy; j++){
			for( dint k = 1; k < mesh->ncellsz; k++){

				dint cell = computeLinearIndex(i,j,k, mesh->ncellsx,mesh->ncellsy,mesh->ncellsz);
				int fid   = mesh->faceACGE + cell*mesh->ncellFaces;

				PetscInt row_id    = cell;
				PetscInt col_id    = mesh->faceToNeighborCell[fid];
				PetscScalar Avalue = mesh->aN[fid];

				ierr = MatSetValues(A,1,&row_id,1,&col_id,&Avalue,INSERT_VALUES);CHKERRQ(ierr);

			}
		}
	}



	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	//	MatView(A, viewer);

	// Source Vector B Assemble Values
	for(int cell = 0;cell < mesh->nicells; ++cell){

		PetscInt row_id = cell;
		PetscScalar Bvalue = mesh->bP[cell];

		ierr = VecSetValues(B,1,&row_id,&Bvalue,INSERT_VALUES);CHKERRQ(ierr);
	}



	ierr = VecAssemblyBegin(B);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(B);CHKERRQ(ierr);

	//    VecView(B,viewer);

	// Creating linear solver and preconditioner
	ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr); // Creates the default KSP context
	ierr = KSPSetType(ksp, KSPBCGS);
	ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr); // Sets the matrix associated with the linear system and a (possibly) different one associated with the preconditioner

	ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
	// ierr = PCSetType(pc,PCNONE);CHKERRQ(ierr);
	ierr = PCSetType(pc,PCJACOBI);CHKERRQ(ierr);
	ierr = KSPSetTolerances(ksp,1.e-20,1.e-20,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);

	ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);  // Sets KSP options from the options database

	ierr = KSPSolve(ksp,B,X);CHKERRQ(ierr);

	// View solver info
   //	ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

	ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"\t Number of Iterations =  %d\n",its);CHKERRQ(ierr);

	// Copy the data to vector
    for(PetscInt i = 0;i < mesh->nicells; ++i){

    	PetscScalar value; // Value of Source Vector

    	ierr = VecGetValues(X,1,&i,&value); CHKERRQ(ierr);

    	int cell = i;

    	alpha[cell] =  value;
    }



//    VecView(X,viewer);

    ierr = MatDestroy(&A);CHKERRQ(ierr);
    ierr = VecDestroy(&X);CHKERRQ(ierr);
    ierr = VecDestroy(&B);CHKERRQ(ierr);
    ierr = KSPDestroy(&ksp);CHKERRQ(ierr);


	return ierr;

}

