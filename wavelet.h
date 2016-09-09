
#ifndef WAVELET_H
#define WAVELET_H


#include "exports.h"

/* Discrete Wavelet Transform
input parameters:
	pImage: image buffer
	ht,wd:  input image height and width
output parameters:
	pDWT: buffer to save wavelet decomposition image
	dstHt, dstWd: width and height of wavelet decomposition image
*/
DLL_EXPORT void DWTProcess(unsigned char* pImage, int ht, int wd, float** pDWT, int* dstHt, int* dstWd, int level);



#endif