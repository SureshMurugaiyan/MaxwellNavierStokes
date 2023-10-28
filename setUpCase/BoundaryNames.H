#ifndef __BOUNDARYNAMES_H
#define __BOUNDARYNAMES_H

typedef enum {
	none				= 0,
	topBoundary			= 1,
	bottomBoundary		= 2,
	leftBoundary		= 3,
	rightBoundary		= 4,
	frontBoundary       = 5,
	backBoundary        = 6,
	inlet				= 7,
	outlet				= 8,
} variable_name;

#endif
