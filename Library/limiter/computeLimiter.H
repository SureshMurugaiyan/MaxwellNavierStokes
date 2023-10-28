#include "mesh.h"

//#include "../Library/limiter/computeLimiter.H"

dfloat computeLimiter
(
	dfloat Flux,
	dfloat phiP,
	dfloat phiN,
	dfloat* gradCP,
	dfloat* gradCN,
	dfloat* d,
	dint interpolationScheme,
	dfloat blendingCoeff,
	mesh_t*	mesh
);

dfloat computeLimiter
(
	dfloat Flux,
	dfloat phiP,
	dfloat phiN,
	dfloat* gradCP,
	dfloat* gradCN,
	dfloat* d,
	dint interpolationScheme,
	dfloat blendingCoeff,
	mesh_t*	mesh
)
{

	blendingCoeff = stabilise(blendingCoeff, mesh); // avoid division by zero

	dfloat lambda  = 0.0;

	if (interpolationScheme == linear) {

		lambda  = 1.0;
	}

	else if(interpolationScheme == upwind) {

		lambda  = 0.0;

	}

	else if (interpolationScheme == interfaceCompression ){

		lambda  = mymin(mymax(1 - mymax(sqr(1 - 4*phiP*(1 - phiP)), sqr(1 - 4*phiN*(1 - phiN))),0), 1);

	}

	else if (interpolationScheme == Gamma ){

		dfloat phiC = computePhiC(Flux,phiP,phiN,gradCP,gradCN,d,mesh);
		lambda 		= mymin(mymax(phiC/blendingCoeff, 0), 1);
	}

	else if (interpolationScheme == VanLeer ){

		dfloat r = compute_r(Flux,phiP,phiN,gradCP,gradCN,d,mesh);
		lambda 	= (r + magnitude(r))/(1 + magnitude(r));
	}

	else if (interpolationScheme == MUSCL ){

		dfloat r = compute_r(Flux,phiP,phiN,gradCP,gradCN,d,mesh);
		lambda 	= mymax(mymin(mymin(2*r, 0.5*r + 0.5), 2), 0);
	}

	else if (interpolationScheme == Minmod ){

		dfloat r = compute_r(Flux,phiP,phiN,gradCP,gradCN,d,mesh);
		lambda 	= mymax(mymin(r, 1), 0);
	}

	else if (interpolationScheme ==  SuperBee){

		dfloat r = compute_r(Flux,phiP,phiN,gradCP,gradCN,d,mesh);
		lambda 	= mymax(mymax(mymin(2*r, 1), mymin(r, 2)), 0);
	}

	return lambda;


}
