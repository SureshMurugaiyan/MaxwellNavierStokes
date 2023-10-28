#include "mesh.h"

//#include "../Library/limiter/computePhiC.H"

dfloat computePhiC
(
	dfloat Flux,
	dfloat phiP,
	dfloat phiN,
	dfloat* gradCP,
	dfloat* gradCN,
	dfloat* d,
	mesh_t*	mesh
);

dfloat computePhiC
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

	gradphiCdotd = stabilise(gradphiCdotd, mesh);

	dfloat phiC = 1.0 - 0.5*gradphifdotd/gradphiCdotd;

	return phiC;


}
