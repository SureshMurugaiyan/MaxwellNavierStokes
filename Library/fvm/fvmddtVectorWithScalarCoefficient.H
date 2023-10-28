#include "mesh.h"

//#include "../Library/fvm/fvmddtVectorWithScalarCoefficient.H"


dfloat fvmddtVectorWithScalarCoefficient
(
	dfloat GT,
	dfloat* alpha,
	dfloat* Phi,
	dint PhiName,
	dint timeScheme,
	mesh_t *mesh
);

dfloat fvmddtVectorWithScalarCoefficient
(
	dfloat GT,
	dfloat* alpha,
	dfloat* Phi,
	dint PhiName,
	dint timeScheme,
	mesh_t *mesh
)
{

	if (mesh->timeStep==0) {
		timeScheme	= euler; // use  euler for first time step ,
	}



	/*************************************************************/
	/* EULER  SCHEME                                             */
	/*************************************************************/

	if (timeScheme	== euler){

		for (dint cell = 0; cell < mesh->ncells; cell++) {

			dint cell_O  = cell + mesh->ncells;
			dint cell_OO = cell + mesh->ncells*2;

			for (dint dimension = 0; dimension < mesh->ndimension; dimension++) {

				dint cell_dim 		= dimension + cell*mesh->ndimension;
				dint cell_O_dim 	= dimension + cell_O*mesh->ndimension;
				dint cell_OO_dim 	= dimension + cell_OO*mesh->ndimension;

				dfloat AP_new  = GT*alpha[cell]/mesh->deltatime;
				dfloat BP_new  = GT*alpha[cell]*Phi[cell_O_dim]/mesh->deltatime;

				mesh->AP[cell_dim]  =  mesh->AP[cell_dim] + AP_new;
				mesh->BP[cell_dim]  =  mesh->BP[cell_dim] + BP_new;


			}

		}

	}


	/*************************************************************/
	/* EULER BACKWARD SCHEME                                     */
	/*************************************************************/

	if (timeScheme	== eulerBackward){

		for (dint cell = 0; cell < mesh->ncells; cell++) {

			dint cell_O  = cell + mesh->ncells;
			dint cell_OO = cell + mesh->ncells*2;

			for (dint dimension = 0; dimension < mesh->ndimension; dimension++) {

				dint cell_dim 		= dimension + cell*mesh->ndimension;
				dint cell_O_dim 	= dimension + cell_O*mesh->ndimension;
				dint cell_OO_dim 	= dimension + cell_OO*mesh->ndimension;

				dfloat AP_new  = GT*alpha[cell]*3.0/(2.0*mesh->deltatime);
				dfloat BP_new  = GT*alpha[cell]*(4.0*Phi[cell_O_dim] - Phi[cell_OO_dim])/(2.0*mesh->deltatime);

				mesh->AP[cell_dim]  =  mesh->AP[cell_dim] + AP_new;
				mesh->BP[cell_dim]  =  mesh->BP[cell_dim] + BP_new;


			}

		}

	}



	return 0.0;

}
