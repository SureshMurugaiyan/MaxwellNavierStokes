#include "mesh.h"

void createSurfaceMesh(mesh_t *mesh);

void createSurfaceMesh(mesh_t *mesh){

	createBottomFaceSurfaceMesh(mesh);

	createTopFaceSurfaceMesh(mesh);

	createBackFaceSurfaceMesh(mesh);

	createFrontFaceSurfaceMesh(mesh);

	createLeftFaceSurfaceMesh(mesh);

	createRightFaceSurfaceMesh(mesh);

//	printf ( "!---------------------------------------------------------------------------------------------------- !\n" );
//	printf ( "! generate Surface Mesh                                                                               !\n" );
//	printf ( "!---------------------------------------------------------------------------------------------------- !\n" );
//
//	printVertexVector(mesh->vertices,mesh);

}
