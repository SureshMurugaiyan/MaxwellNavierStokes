#include "mesh.h"

void setModifiedInnerOuterRadius(mesh_t *mesh);

void setModifiedInnerOuterRadius(mesh_t *mesh){

dfloat	Hi = mesh->innerRadius*tan(mesh->axisSymmetricMeshAngle/2.0);
dfloat	Ho = mesh->outerRadius*tan(mesh->axisSymmetricMeshAngle/2.0);

mesh->innerRadius = sqrt(mesh->innerRadius*mesh->innerRadius+ Hi*Hi);
mesh->outerRadius = sqrt(mesh->outerRadius*mesh->outerRadius+ Ho*Ho);

}


