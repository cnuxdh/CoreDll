

#ifndef MOSAICFUNCS_H
#define MOSAICFUNCS_H

#include "exports.h"


#include "cv.h"
#include "highgui.h"
#include "cxcore.h"

#include "../../Corelib/commondata.h"

int FundamentalRefine(MyPointF* pLeftPt, MyPointF* pRightPt, int npt, CvMat** fundamental_matrix);

int DeleteRedundantPairs(MyPointF* pt1, MyPointF* pt2, int np);


int DeleteSamePoints(struct feature* feat, int n);


void LoadSiftFile(char* filename, struct feature** feat, int* n);

IplImage* stack_imgs( IplImage* img1, IplImage* img2 );

IplImage* MosaicTwoImagesPolynomial(IplImage* srcImage, IplImage* dstImage, float* p);

IplImage* MosaicTwoImages(IplImage* srcImage, IplImage* dstImage, float* p, int* minx, int* miny);


void DLL_EXPORT FindCorrespondence(struct feature* feat1, int n1,
						struct feature* feat2, int n2,
						struct feature* p1, struct feature* p2, int* np);

void DLL_EXPORT FindCorrespondence(struct feature* feat1, int n1,
						struct feature* feat2, int n2,
						MyPointF* p1, MyPointF* p2, int* np);


void FindCorrespondenceNew(struct feature* feat1, int n1,
						   struct feature* feat2, int n2,
						   MyPointF* p1, MyPointF* p2, 
						   int *ip1, int* ip2,
						   int* np);

void GenerateTriangle(IplImage* img, MyPointF* pt, int np);


IplImage*  MosaicUsingTriangle(IplImage* srcImage, IplImage* dstImage, 
						 MyPointF* srtPt, MyPointF* dstPt, int np);


IplImage*  MosaicUsingHomography(IplImage* srcImage, IplImage* dstImage, double* H, int* nminx=NULL, int* miny=NULL);

//特征点加密
//加密后的点存放到原来的数组中
void  AddFeaturePoint(IplImage* srcImage, IplImage* dstImage, MyPointF** srcPt, MyPointF** dstPt, int* np);


//mosaic based on bundle result
int DLL_EXPORT BundleMosaic(char* listfile, char* bundleFile, int ht, int wd, const char* savePath=NULL);
int DLL_EXPORT BundlerMosaicOnebyOne(char* listfile, char* bundleFile, int ht, int wd, const char* savePath=NULL);

int  DLL_EXPORT  BundlerMosaicAR(char* listfile, char* bundleFile, int ht, int wd, const char* savePath=NULL);
void DLL_EXPORT  FitPlane1(double* px, double* py, double* pz, int np, double* R);

int PlaneMosaic(char* listfile, char* bundleFile, double ratio, const char* savePath);



#endif