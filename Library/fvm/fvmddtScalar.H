#include "mesh.h"

//#include "../Library/fvm/fvmddtScalar.H"

dfloat fvmDdtOfScalar
(
	dfloat GT,
	dfloat* Phi,
	dint PhiName,
	dint timeScheme,
	mesh_t *mesh
);

dfloat fvmDdtOfScalar
(
	dfloat GT,
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

		for (dint cell = 0; cell < mesh->nicells; cell++) {

			dint cell_O  = cell + mesh->ncells;
			dint cell_OO = cell + mesh->ncells*2;

			dfloat AP_new = GT/mesh->deltatime;
			dfloat R_new  = GT*Phi[cell_O]/mesh->deltatime;

			mesh->aP[cell] = mesh->aP[cell] + AP_new;
			mesh->bP[cell] = mesh->bP[cell] + R_new;

		}

	}

	/*************************************************************/
	/* EULER BACKWARD SCHEME                                     */
	/*************************************************************/

	if (timeScheme	== eulerBackward){

		for (dint cell = 0; cell < mesh->nicells; cell++) {

			dint cell_O  = cell + mesh->ncells;
			dint cell_OO = cell + mesh->ncells*2;

			dfloat AP_new = GT*3.0/(2.0*mesh->deltatime);
			dfloat R_new  = GT*(4.0*Phi[cell_O] - Phi[cell_OO])/(2.0*mesh->deltatime);

			mesh->aP[cell] = mesh->aP[cell]  + AP_new;
			mesh->bP[cell] = mesh->bP[cell]  + R_new;

		}
	}


	return 0.0;

}
