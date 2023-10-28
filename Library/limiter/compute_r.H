#include "mesh.h"

//#include "../Library/limiter/compute_r.H"

dfloat compute_r // Compute r ratio
(
	dfloat Flux,
	dfloat phiP,
	dfloat phiN,
	dfloat* gradCP,
	dfloat* gradCN,
	dfloat* d,
	mesh_t*	mesh
);

dfloat compute_r // Compute r ratio
(
	dfloat Flux,
	dfloat phiP,
	dfloat phiN,
	dfloat* gradCP,
	dfloat* gradCN,
	dfloat* d,
	mesh_t*	mesh
)
{

	dfloat gradphifdotd = 0.0;
	dfloat gradphiCdotd = 0.0;


	gradphifdotd = phiN - phiP;


	if(Flux > 0){

		gradphiCdotd = dotProduct(d,gradCP);
	}
	else{
		gradphiCdotd = dotProduct(d,gradCN);
	}

	gradphifdotd = stabilise(gradphifdotd, mesh);

	dfloat r = 2.0*gradphiCdotd/gradphifdotd-1.0;

	return r;


}
