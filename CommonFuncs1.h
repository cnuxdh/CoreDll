#ifndef COMMON_FUNCS
#define COMMON_FUNCS

#include "corelib/commondata.h"
//#include "exports.h"

#include "cv.h"
#include "highgui.h"

#include <vector>
using namespace std;


#define DLL_EXPORT  _declspec(dllexport)


double DLL_EXPORT BilinearPixel(double ix, double iy, IplImage* pImage);

void   DLL_EXPORT IplImageToGrayImage(IplImage* pImage, unsigned char** pBuffer, int* ht, int* wd);
void   DLL_EXPORT IplImageToColorImage(IplImage* pImage, unsigned char** r, unsigned char** g, unsigned char** b, int* ht, int* wd);
void   DLL_EXPORT IplImageToFloatImage(IplImage* pImage, float** pBuffer, int* ht, int* wd);

void   DLL_EXPORT GenerateOpticalFlow(IplImage* pFirst, IplImage* pSecond, char* filepath);

void   DLL_EXPORT ConnectedComonent(IplImage* pImage, MyRect* pRect, int* nRect);

void   DLL_EXPORT VerticalHaarSeg(IplImage* pSrc, IplImage* pDst);

void   DLL_EXPORT PlateImageSeg(IplImage* pImage, unsigned char* mask);

void   DLL_EXPORT VerticalEdgeDetect(IplImage* pSrc, IplImage* pDst);

void   DLL_EXPORT IplImageSplit(IplImage* pSrc, IplImage* pR, IplImage* pG, IplImage* pB);

void   DLL_EXPORT MergeParallelLines(vector<stLINE>& lines);



#endif