#include "mesh.h"

dint computeLocalFace(dint fid, mesh_t *mesh);

dint computeLocalFace(dint fid, mesh_t *mesh){

	dint LocalFace = fid%mesh->ncellFaces;

	return  LocalFace;


}


