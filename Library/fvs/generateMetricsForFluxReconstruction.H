#include "mesh.h"

//#include "../Library/fvs/generateMetricsForFluxReconstruction.H"

void generateMetricsForFluxReconstruction(mesh_t *mesh);

void generateMetricsForFluxReconstruction(mesh_t *mesh){

	// Compute SnDyadSf at cell Faces

	for (dint cell = 0; cell < mesh->nicells; ++cell) {

		for (dint face = 0; face < mesh->ncellFaces; ++face) {

			dint fid	=  face  + cell*mesh->ncellFaces;

			for(dint i = 0; i < mesh->ndimension; i++){

				for(dint j = 0; j < mesh->ndimension; j++){

					dint i_dim  = i + fid*mesh->ndimension;
					dint j_dim  = j + fid*mesh->ndimension;

					dint id = j + i*mesh->ndimension + fid*mesh->ndimension*mesh->ndimension;

					mesh->SnDyadSf[id] = mesh->Sn[i_dim]*mesh->Sf[j_dim];


				}

			}

		}

	}

	//  Compute Summed SnDyadSf  cell Faces gives cell center tensor

	dfloat *SnDyadSf;
	SnDyadSf = (dfloat*) calloc(mesh->ndimension*mesh->ndimension, sizeof(dfloat));


	for (dint cell = 0; cell < mesh->nicells; ++cell) {

		for (dint n = 0; n < mesh->ndimension*mesh->ndimension; ++n) {SnDyadSf[n]=0;}


		for (dint face = 0; face < mesh->ncellFaces; ++face) {

			dint fid	=  face  + cell*mesh->ncellFaces;

			for(dint i = 0; i < mesh->ndimension; i++){

				for(dint j = 0; j < mesh->ndimension; j++){

					dint i_dim  = i + fid*mesh->ndimension;
					dint j_dim  = j + fid*mesh->ndimension;

					dint id = j + i*mesh->ndimension + fid*mesh->ndimension*mesh->ndimension;

					dint idt = j+i*mesh->ndimension;

					SnDyadSf[idt] = SnDyadSf[idt] + mesh->SnDyadSf[id];

				}

			}

		}


		// Store the summed face values at cell centers

		for(dint i = 0; i < mesh->ndimension; i++){

			for(dint j = 0; j < mesh->ndimension; j++){

				dint id = j + i*mesh->ndimension + cell*mesh->ndimension*mesh->ndimension;

				dint idt = j+i*mesh->ndimension;

				mesh->summedSnDyadSf[id] = SnDyadSf[idt] ;

			}

		}


	}

	// Computing inverse of cell center tensor field

	dint m,n;
	dint a,b,c,d,e,f,g,h,i;
	dint A,B,C,D,E,F,G,H,I;

	for (dint cell = 0; cell < mesh->nicells; ++cell) {

		m=0;n=0;
		a = n + m*mesh->ndimension + cell*mesh->ndimension*mesh->ndimension;

		m=0;n=1;
		b = n + m*mesh->ndimension + cell*mesh->ndimension*mesh->ndimension;

		m=0;n=2;
		c = n + m*mesh->ndimension + cell*mesh->ndimension*mesh->ndimension;

		m=1;n=0;
		d = n + m*mesh->ndimension + cell*mesh->ndimension*mesh->ndimension;

		m=1;n=1;
		e = n + m*mesh->ndimension + cell*mesh->ndimension*mesh->ndimension;

		m=1;n=2;
		f = n + m*mesh->ndimension + cell*mesh->ndimension*mesh->ndimension;

		m=2;n=0;
		g = n + m*mesh->ndimension + cell*mesh->ndimension*mesh->ndimension;

		m=2;n=1;
		h = n + m*mesh->ndimension + cell*mesh->ndimension*mesh->ndimension;

		m=2;n=2;
		i = n + m*mesh->ndimension + cell*mesh->ndimension*mesh->ndimension;

		//====================!

		A = a;
		B = d;
		C = g;
		D = b;
		E = e;
		F = h;
		G = c;
		H = f;
		I = i;

		dfloat aa = mesh->summedSnDyadSf[a];
		dfloat bb = mesh->summedSnDyadSf[b];
		dfloat cc = mesh->summedSnDyadSf[c];
		dfloat dd = mesh->summedSnDyadSf[d];
		dfloat ee = mesh->summedSnDyadSf[e];
		dfloat ff = mesh->summedSnDyadSf[f];
		dfloat gg = mesh->summedSnDyadSf[g];
		dfloat hh = mesh->summedSnDyadSf[h];
		dfloat ii = mesh->summedSnDyadSf[i];


		dfloat AA =  (ee*ii - ff*hh);
		dfloat BB = -(dd*ii - ff*gg);
		dfloat CC =  (dd*hh - ee*gg);
		dfloat DD = -(bb*ii - cc*hh);
		dfloat EE =  (aa*ii - cc*gg);
		dfloat FF = -(aa*hh - bb*gg);
		dfloat GG =  (bb*ff - cc*ee);
		dfloat HH = -(aa*ff - cc*dd);
		dfloat II =  (aa*ee - bb*dd);

		dfloat det = aa*AA + bb*BB + cc*CC;
		dfloat rdet = 1.0/stabilise(det,mesh);

		mesh->inverseOfSummedSnDyadSf[A] = AA*rdet ;
		mesh->inverseOfSummedSnDyadSf[B] = BB*rdet ;
		mesh->inverseOfSummedSnDyadSf[C] = CC*rdet ;
		mesh->inverseOfSummedSnDyadSf[D] = DD*rdet ;
		mesh->inverseOfSummedSnDyadSf[E] = EE*rdet ;
		mesh->inverseOfSummedSnDyadSf[F] = FF*rdet ;
		mesh->inverseOfSummedSnDyadSf[G] = GG*rdet ;
		mesh->inverseOfSummedSnDyadSf[H] = HH*rdet ;
		mesh->inverseOfSummedSnDyadSf[I] = II*rdet ;


//		// redo for 2 by 2 matrix for 2D since 3D causes determinant zero
//
//		det = aa*ee - bb*dd;
//		rdet = 1.0/stabilise(det,mesh);
//
//		mesh->inverseOfSummedSnDyadSf[A] =  ee*rdet ;
//		mesh->inverseOfSummedSnDyadSf[B] = -dd*rdet ;
//		mesh->inverseOfSummedSnDyadSf[C] =  0 ;
//		mesh->inverseOfSummedSnDyadSf[D] = -bb*rdet ;
//		mesh->inverseOfSummedSnDyadSf[E] =  aa*rdet ;
//		mesh->inverseOfSummedSnDyadSf[F] =  0 ;
//		mesh->inverseOfSummedSnDyadSf[G] =  0 ;
//		mesh->inverseOfSummedSnDyadSf[H] =  0 ;
//		mesh->inverseOfSummedSnDyadSf[I] =  0 ;


//		if( cell ==0){
//			printf ( " (%9.2e) ", mesh->inverseOfSummedSnDyadSf[A] );
//			printf ( " (%9.2e) ", mesh->inverseOfSummedSnDyadSf[D] );
//			printf ( " (%9.2e) \n", mesh->inverseOfSummedSnDyadSf[G] );
//			printf ( " (%9.2e) ", mesh->inverseOfSummedSnDyadSf[B] );
//			printf ( " (%9.2e) ", mesh->inverseOfSummedSnDyadSf[E] );
//			printf ( " (%9.2e) \n", mesh->inverseOfSummedSnDyadSf[H] );
//			printf ( " (%9.2e) ", mesh->inverseOfSummedSnDyadSf[C] );
//			printf ( " (%9.2e) ", mesh->inverseOfSummedSnDyadSf[F] );
//			printf ( " (%9.2e) \n", mesh->inverseOfSummedSnDyadSf[I] );
//		}

//		if( cell ==0){
//			printf ( " (%9.2e) ", mesh->summedSnDyadSf[a] );
//			printf ( " (%9.2e) ", mesh->summedSnDyadSf[b] );
//			printf ( " (%9.2e) \n", mesh->summedSnDyadSf[c] );
//			printf ( " (%9.2e) ", mesh->summedSnDyadSf[d] );
//			printf ( " (%9.2e) ", mesh->summedSnDyadSf[e] );
//			printf ( " (%9.2e) \n", mesh->summedSnDyadSf[f] );
//			printf ( " (%9.2e) ", mesh->summedSnDyadSf[g] );
//			printf ( " (%9.2e) ", mesh->summedSnDyadSf[h] );
//			printf ( " (%9.2e) \n", mesh->summedSnDyadSf[i] );
//		}




	}




	free(SnDyadSf);

}
