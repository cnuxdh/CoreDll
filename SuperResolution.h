#ifndef SUPER_RESOLUTION_H
#define SUPER_RESOLUTION_H

#include "exports.h"


/*
   x,y: aligned point ordinates
   v:   pixel color value
   npt: number of points
*/
DLL_EXPORT void POCS(double* x, double* y, int* v, int npt, char* outfile);





#endif