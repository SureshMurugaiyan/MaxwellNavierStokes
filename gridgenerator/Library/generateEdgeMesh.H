#include "mesh.h"

void generateEdgeMesh(mesh_t *mesh);

void generateEdgeMesh(mesh_t *mesh){

	// Corner Edges

	generateNorthWestCornerEdge(mesh);

	generateSouthWestCornerEdge(mesh);

	generateSouthEastCornerEdge(mesh);

	generateNorthEastCornerEdge(mesh);


	// Bottom Edges

	generateBottomNorthEdgeMesh(mesh);

	generateBottomSouthEdgeMesh(mesh);

	generateBottomWestEdgeMesh(mesh);

	generateBottomEastEdgeMesh(mesh);


	// Top Edges

	generateTopNorthEdgeMesh(mesh);

	generateTopSouthEdgeMesh(mesh);

	generateTopWestEdgeMesh(mesh);

	generateTopEastEdgeMesh(mesh);


//	printf ( "!--------------------------------------------- !\n" );
//	printf ( "! generate Edge Mesh                           !\n" );
//	printf ( "!--------------------------------------------- !\n" );
//
//	printVertexVector(mesh->vertices,mesh);


}
