#include "mesh.h"

void createTopFaceSurfaceMesh(mesh_t *mesh);

void createTopFaceSurfaceMesh(mesh_t *mesh){


	// Top Face (X8,Y8,Z8),(X5,Y5,Z5),(X6,Y6,Z6) and (X7,Y7,Z7)
	// i = 0                    --> West Boundary
	// i = (mesh->nverticesx-1) --> East Boundary
	// j = 0                    --> South Boundary
	// j = (mesh->nverticesy-1) --> North Boundary
	// k = 0                    --> Bottom Plane
	// k = (mesh->nverticesz-1) --> Top Plane


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

	// X(xi), Y(eta), Z(zeta)

	linearIndex4D_t *linearIndex4D = (linearIndex4D_t*) calloc(1, sizeof(mesh_t));


	dint k = mesh->nverticesz-1; // Top Face

	for (dint i = 1;i < (mesh->nverticesx-1);i++){ // Skipping the first and last corner points
		for (dint j = 1;j < (mesh->nverticesy-1);j++){ // Skipping the first and last corner points
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

				linearIndex4D->i = 0;
				linearIndex4D->j = (mesh->nverticesy-1);
				linearIndex4D->k = k;
				linearIndex4D->l = dim;
				dint nw_vid_dim  = computeLinearIndex4D(linearIndex4D);

				linearIndex4D->i = 0;
				linearIndex4D->j = 0;
				linearIndex4D->k = k;
				linearIndex4D->l = dim;
				dint sw_vid_dim  = computeLinearIndex4D(linearIndex4D);

				linearIndex4D->i = (mesh->nverticesx-1);
				linearIndex4D->j = 0;
				linearIndex4D->k = k;
				linearIndex4D->l = dim;
				dint se_vid_dim  = computeLinearIndex4D(linearIndex4D);

				linearIndex4D->i = (mesh->nverticesx-1);
				linearIndex4D->j = (mesh->nverticesy-1);
				linearIndex4D->k = k;
				linearIndex4D->l = dim;
				dint ne_vid_dim  = computeLinearIndex4D(linearIndex4D);

				linearIndex4D->i = i;
				linearIndex4D->j = (mesh->nverticesy-1);
				linearIndex4D->k = k;
				linearIndex4D->l = dim;
				dint n_vid_dim   = computeLinearIndex4D(linearIndex4D);

				linearIndex4D->i = i;
				linearIndex4D->j = 0;
				linearIndex4D->k = k;
				linearIndex4D->l = dim;
				dint s_vid_dim   = computeLinearIndex4D(linearIndex4D);

				linearIndex4D->i = 0;
				linearIndex4D->j = j;
				linearIndex4D->k = k;
				linearIndex4D->l = dim;
				dint w_vid_dim   = computeLinearIndex4D(linearIndex4D);

				linearIndex4D->i = (mesh->nverticesx-1);
				linearIndex4D->j = j;
				linearIndex4D->k = k;
				linearIndex4D->l = dim;
				dint e_vid_dim   = computeLinearIndex4D(linearIndex4D);

				mesh->vertices[vid_dim] = (1-eta)*mesh->vertices[s_vid_dim]      + eta*mesh->vertices[n_vid_dim] +
						                  (1- xi)*mesh->vertices[w_vid_dim]      +  xi*mesh->vertices[e_vid_dim] -
										  (xi*eta*mesh->vertices[ne_vid_dim]     +  xi*(1-eta)*mesh->vertices[se_vid_dim] +
										   eta*(1-xi)*mesh->vertices[nw_vid_dim] + (1-xi)*(1-eta)*mesh->vertices[sw_vid_dim]);

			}
		}
	}



	free(linearIndex4D);

}
