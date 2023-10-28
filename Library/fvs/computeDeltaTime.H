#include "mesh.h"

//#include "../Library/fvs/computeDeltaTime.H"

void computeDeltaTime(mesh_t *mesh);

void computeDeltaTime(mesh_t *mesh){

	computeCourantNumberUsingConvectiveVelocity(mesh);


	mesh->Co0 = mesh->CoU;


	if(mesh->timeStepSizeControl==targetCourantNumber){

		mesh->deltatt = mesh->Cot*mesh->deltatime/stabilise(mesh->Co0,mesh);

		if(mesh->Co0 >(mesh->Kl*mesh->Cot)){
			mesh->deltatime =mesh->deltatt ;
		}else{
			mesh->deltatime = mymin(mesh->deltatime+mesh->lambdat*(mesh->deltatt-mesh->deltatime) ,2*mesh->deltatime) ;
		}

	}

}
