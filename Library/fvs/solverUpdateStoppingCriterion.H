#include "mesh.h"

//#include "../Library/fvs/solverUpdateStoppingCriterion.H"

void solverUpdateStoppingCriterion(mesh_t *mesh);

void solverUpdateStoppingCriterion(mesh_t *mesh){

	  if (mesh->timeStep == mesh->timeStepMax) {
		  mesh->stopCriterion = true;
	  }

	  if (mesh->time > mesh->finalTime) {
		  mesh->stopCriterion = true;
	  }

	  if (mesh->timeStepMax == 0) {
		  mesh->stopCriterion = true;
	  }

}
