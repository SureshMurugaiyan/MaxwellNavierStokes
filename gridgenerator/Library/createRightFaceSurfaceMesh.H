#include "mesh.h"

void createRightFaceSurfaceMesh(mesh_t *mesh);

void createRightFaceSurfaceMesh(mesh_t *mesh){


	// Right Face (X2,Y2,Z2), (X3,Y3,Z3), (X7,Y7,Z7) and (X6,Y6,Z6)

	// i = 0                    --> West Boundary
	// i = (mesh->nverticesx-1) --> East Boundary
	// j = 0                    --> South Boundary
	// j = (mesh->nverticesy-1) --> North Boundary
	// k = 0                    --> Bottom Plane
	// k = (mesh->nverticesz-1) --> Top Plane

	// BOTTOM BOUNDARY ( Z axis is out of Plane )
	/*        	 Y(j)
	 *  	  	 ^
	 *  	  	 |
	 * 	         |			NORTH
	 * (X4,Y4,Z4)*--------------------------*(X3,Y3,Z3)
	 *  	     |                          |
	 * 	     W	 |                          |  E
	 *       E   |			                |  A
	 *       S   |                          |  S
	 *       T	 |                          |  T
	 *  	     |                          |
	 * (X1,Y1,Z1)*--------------------------*(X2,Y2,Z2)-----------> X(i)
	 * 				       SOUTH
	 */

	// TOP BOUNDARY ( Z axis is out of Plane )
	/*        	 Y(j)
	 *  	  	 ^
	 *  	  	 |
	 * 	         |			NORTH
	 * (X8,Y8,Z8)*--------------------------*(X7,Y7,Z7)
	 *  	     |                          |
	 * 	     W	 |                          |  E
	 *       E   |			                |  A
	 *       S   |                          |  S
	 *       T	 |                          |  T
	 *  	     |                          |
	 * (X5,Y5,Z5)*--------------------------*(X6,Y6,Z6)-----------> X(i)
	 * 				       SOUTH
	 */

	linearIndex4D_t *linearIndex4D = (linearIndex4D_t*) calloc(1, sizeof(mesh_t));

	dint i = mesh->nverticesx-1;

	for (dint j = 1;j < (mesh->nverticesy-1);j++){ // Skipping the first and last corner points
		for (dint k = 1;k < (mesh->nverticesz-1);k++){ // Skipping the first and last corner points
			for(dint dim = 0; dim < mesh->ndimension; dim++){

				dfloat xi   = mesh->xi[i];
				dfloat eta  = mesh->eta[j];
				dfloat zeta = mesh->zeta[k];

				dint vid     = computeLinearIndex(i,j,k, mesh->nverticesx, mesh->nverticesy, mesh->nverticesz);
				dint vid_dim = dim + vid*mesh->ndimension;

				linearIndex4D->iN = mesh->nverticesx;
				linearIndex4D->jN = mesh->nverticesy;
				linearIndex4D->kN = mesh->nverticesz;
				linearIndex4D->lN = mesh->ndimension;

				linearIndex4D->i = i;
				linearIndex4D->j = 0;
				linearIndex4D->k = (mesh->nverticesz-1);
				linearIndex4D->l = dim;
				dint nw_vid_dim  = computeLinearIndex4D(linearIndex4D);

				linearIndex4D->i = i;
				linearIndex4D->j = 0;
				linearIndex4D->k = 0;
				linearIndex4D->l = dim;
				dint sw_vid_dim  = computeLinearIndex4D(linearIndex4D);

				linearIndex4D->i = i;
				linearIndex4D->j = (mesh->nverticesy-1);
				linearIndex4D->k = 0;
				linearIndex4D->l = dim;
				dint se_vid_dim  = computeLinearIndex4D(linearIndex4D);

				linearIndex4D->i = i;
				linearIndex4D->j = (mesh->nverticesy-1);
				linearIndex4D->k = (mesh->nverticesz-1);
				linearIndex4D->l = dim;
				dint ne_vid_dim  = computeLinearIndex4D(linearIndex4D);

				linearIndex4D->i = i;
				linearIndex4D->j = j;
				linearIndex4D->k = (mesh->nverticesz-1);
				linearIndex4D->l = dim;
				dint n_vid_dim   = computeLinearIndex4D(linearIndex4D);

				linearIndex4D->i = i;
				linearIndex4D->j = j;
				linearIndex4D->k = 0;
				linearIndex4D->l = dim;
				dint s_vid_dim   = computeLinearIndex4D(linearIndex4D);

				linearIndex4D->i = i;
				linearIndex4D->j = 0;
				linearIndex4D->k = k;
				linearIndex4D->l = dim;
				dint w_vid_dim   = computeLinearIndex4D(linearIndex4D);

				linearIndex4D->i = i;
				linearIndex4D->j = (mesh->nverticesy-1);
				linearIndex4D->k = k;
				linearIndex4D->l = dim;
				dint e_vid_dim   = computeLinearIndex4D(linearIndex4D);

				mesh->vertices[vid_dim] =  (1-zeta)*mesh->vertices[s_vid_dim]      + zeta*mesh->vertices[n_vid_dim] +
						                   (1- eta)*mesh->vertices[w_vid_dim]      +  eta*mesh->vertices[e_vid_dim] -
										  (eta*zeta*mesh->vertices[ne_vid_dim]     +  eta*(1-zeta)*mesh->vertices[se_vid_dim] +
										   zeta*(1-eta)*mesh->vertices[nw_vid_dim] + (1-eta)*(1-zeta)*mesh->vertices[sw_vid_dim]);

			}
		}
	}



	free(linearIndex4D);



}
