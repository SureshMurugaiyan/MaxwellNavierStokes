#include "mesh.h"

dfloat magnitude(dfloat a);

dfloat magnitude(dfloat a){

	dfloat c = 0;

	if(a<0)
	{
		c =  -a;
	}
	else
	{
		c = a;
	}

	return c;
}
