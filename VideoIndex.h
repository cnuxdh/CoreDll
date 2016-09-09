#ifndef VIDEO_INDEX_H
#define VIDEO_INDEX_H

#include "exports.h"


#include "cv.h"
#include "highgui.h"
#include "cvtypes.h"

#include <vector>

using namespace std;


/* number of bins of HSV in histogram */
#define NH 10
#define NS 10
#define NV 10

/* max HSV values */
#define H_MAX 360.0
#define S_MAX 1.0
#define V_MAX 1.0

/* low thresholds on saturation and value for histogramming */
#define S_THRESH 0.1
#define V_THRESH 0.2

/* distribution parameter */
#define LAMBDA 30

typedef struct histogram 
{
	//float histo[NH*NS + NV];   /**< histogram array */
	float histo[256];   /**< histogram array */
	int n;                     /**< length of histogram array */
} histogram;


 DLL_EXPORT void  CalculateGrayHistogram(IplImage* gray, histogram* pHist);
 DLL_EXPORT void  normalize_histogram( histogram* histo );
 DLL_EXPORT float histo_dist_sq( histogram* h1, histogram* h2 );
 DLL_EXPORT void  calc_histogram( IplImage* img, histogram* histo );
 DLL_EXPORT IplImage*  bgr2hsv( IplImage* bgr );
 //DLL_EXPORT int   histo_bin( float h, float s, float v );
 DLL_EXPORT int  VideoSearch(IplImage* pFrame, histogram* pDstHist, CvRect* pRect);
 DLL_EXPORT void VideoSearch(IplImage* pFrame, histogram* pDstHist, vector<CvRect>& vecRect );
 DLL_EXPORT void VideoSearchUsingIntegralHist(IplImage* pFrame, histogram* pDstHist, vector<CvRect>& vecRect );



#endif