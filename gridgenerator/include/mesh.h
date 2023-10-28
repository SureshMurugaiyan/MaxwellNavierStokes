#ifndef __MESH_H
#define __MESH_H



#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "math.h"

#define dfloat double
#define dint int

#define mymax(a,b) ( ((a)>(b)) ? (a):(b) )
#define mymin(a,b) ( ((a)<(b)) ? (a):(b) )
#define PI acos(-1)



typedef enum {
	vertexScalar		=	0,
	vertexVector		=	1,
} field_type;

typedef struct{

	//------------------------------------------------------------------!
	// Declarations for setMeshControls
	//------------------------------------------------------------------!

	dint ncellsx; /* number of cells along x-direction*/
	dint ncellsy; /* number of cells along y-direction*/
	dint ncellsz; /* number of cells along z-direction*/

	bool convertCylindricalPolarToCartesian;
	bool convertSphericalPolarToCartesian;

	//------------------------------------------------------------------!
	// Declarations for assignMeshDimensions
	//------------------------------------------------------------------!

	dint ncells;  /* number of cells in mesh */

	dint nverticesx; /* number of vertices along x-direction*/
	dint nverticesy; /* number of vertices along y-direction*/
	dint nverticesz; /* number of vertices along z-direction*/

	dint nvertices;  /* number of vertices in mesh */

	dint ndimension; /* nmeshdimension = x,y,z */

	dint ncellVertices; /*number of vertices per cell */

	//------------------------------------------------------------------!
	// Declarations for calculateIndices
	//------------------------------------------------------------------!

	dint dim_x;
	dint dim_y;
	dint dim_z;

	dint northwesttopVertexID;
	dint southwesttopVertexID;
	dint southeasttopVertexID;
	dint northeasttopVertexID;

	dint northwestbottomVertexID;
	dint southwestbottomVertexID;
	dint southeastbottomVertexID;
	dint northeastbottomVertexID;


	dint northwesttopVertexID_x;
	dint southwesttopVertexID_x;
	dint southeasttopVertexID_x;
	dint northeasttopVertexID_x;

	dint northwesttopVertexID_y;
	dint southwesttopVertexID_y;
	dint southeasttopVertexID_y;
	dint northeasttopVertexID_y;

	dint northwesttopVertexID_z;
	dint southwesttopVertexID_z;
	dint southeasttopVertexID_z;
	dint northeasttopVertexID_z;

	dint northwestbottomVertexID_x;
	dint southwestbottomVertexID_x;
	dint southeastbottomVertexID_x;
	dint northeastbottomVertexID_x;

	dint northwestbottomVertexID_y;
	dint southwestbottomVertexID_y;
	dint southeastbottomVertexID_y;
	dint northeastbottomVertexID_y;

	dint northwestbottomVertexID_z;
	dint southwestbottomVertexID_z;
	dint southeastbottomVertexID_z;
	dint northeastbottomVertexID_z;

	//------------------------------------------------------------------!
	// Declarations for generateComputationalDomain
	//------------------------------------------------------------------!

	dfloat *xi;
	dfloat *eta;
	dfloat *zeta;

	//------------------------------------------------------------------!
	// Declarations
	//------------------------------------------------------------------!

	dfloat *vertices; /* coordinates of vertices in cartesian coordinates*/

	//------------------------------------------------------------------!
	// Declarations for calculateIndices
	//------------------------------------------------------------------!
	dfloat x1;
	dfloat x2;
	dfloat x3;
	dfloat x4;
	dfloat x5;
	dfloat x6;
	dfloat x7;
	dfloat x8;

	dfloat y1;
	dfloat y2;
	dfloat y3;
	dfloat y4;
	dfloat y5;
	dfloat y6;
	dfloat y7;
	dfloat y8;

	dfloat z1;
	dfloat z2;
	dfloat z3;
	dfloat z4;
	dfloat z5;
	dfloat z6;
	dfloat z7;
	dfloat z8;


	//------------------------------------------------------------------!
	// Declarations for writeFiniteVolumeMesh
	//------------------------------------------------------------------!

	dint *cellToVertex; /* Global to Global Vertex connectivity */

	//------------------------------------------------------------------!
	// Declarations for axisymmetric mesh
	//------------------------------------------------------------------!

	dfloat axisSymmetricMeshAngle;

	//------------------------------------------------------------------!
	// Declarations for setUpCase
	//------------------------------------------------------------------!

	dfloat x_min;
	dfloat x_max;

	dfloat y_min;
	dfloat y_max;

	dfloat z_min;
	dfloat z_max;


}mesh_t;

typedef struct{
	dint i;
	dint j;
	dint k;
	dint iN;
	dint jN;
	dint kN;
}linearIndex3D_t;

typedef struct{
	dint i;
	dint j;
	dint k;
	dint l;
	dint iN;
	dint jN;
	dint kN;
	dint lN;
}linearIndex4D_t;


void setMeshControls(mesh_t *mesh);

void assignMeshDimensions(mesh_t *mesh);

void calculateIndices(mesh_t *mesh);

dint computeLinearIndex(dint i, dint j, dint k, dint nx, dint ny, dint nz);

void createField(dfloat** phi, dint variableType, mesh_t *mesh);

void setMeshCorners(mesh_t *mesh);

void assignMeshCorners(mesh_t *mesh);

void generateComputationalDomain(mesh_t *mesh);

void printVertexVector(dfloat* Phi, mesh_t *mesh);

dint computeLinearIndex4D(linearIndex4D_t *linearIndex4D);

void convertCylindricalPolarToCartesian(mesh_t *mesh);

void convertSphericalPolarToCartesian(mesh_t *mesh);

void writeFiniteDifferenceMesh(mesh_t *mesh);

void writeFiniteVolumeMesh(mesh_t *mesh);

void meshCellToVerticesConnectivity(mesh_t *mesh);

void meshConstructComputationalDomainFields(mesh_t *mesh);

void meshDestructComputationalDomainFields(mesh_t *mesh);

void generateNorthWestCornerEdge(mesh_t *mesh);
void generateSouthWestCornerEdge(mesh_t *mesh);
void generateSouthEastCornerEdge(mesh_t *mesh);
void generateNorthEastCornerEdge(mesh_t *mesh);

void generateBottomNorthEdgeMesh(mesh_t *mesh);
void generateBottomSouthEdgeMesh(mesh_t *mesh);
void generateBottomWestEdgeMesh(mesh_t *mesh);
void generateBottomEastEdgeMesh(mesh_t *mesh);

void generateTopNorthEdgeMesh(mesh_t *mesh);
void generateTopSouthEdgeMesh(mesh_t *mesh);
void generateTopWestEdgeMesh(mesh_t *mesh);
void generateTopEastEdgeMesh(mesh_t *mesh);

void generateEdgeMesh(mesh_t *mesh);

void createBottomFaceSurfaceMesh(mesh_t *mesh);
void createTopFaceSurfaceMesh(mesh_t *mesh);
void createBackFaceSurfaceMesh(mesh_t *mesh);
void createFrontFaceSurfaceMesh(mesh_t *mesh);
void createLeftFaceSurfaceMesh(mesh_t *mesh);
void createRightFaceSurfaceMesh(mesh_t *mesh);

void createSurfaceMesh(mesh_t *mesh);

void createVolumeMesh(mesh_t *mesh);


#endif
