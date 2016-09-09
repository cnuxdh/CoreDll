/*
  this file contains the general program definitions.
  revised of Rob Hess's version 1.0.0-20060306
*/

#ifndef DEFS_H
#define DEFS_H

#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
//#include <iostream.h>

/* OpenCV library */
//#include "cv.h"
//#include "highgui.h"
//#include "cxcore.h"


/******************************* Defs and macros *****************************/
#define SIFT_INTVLS				5		// intervals per octave 
#define SIFT_SIGMA				1.6		// sigma for gaussian smoothing 
#define SIFT_CONTR_THR			0.04	// threshold on keypoint contrast |D(x)|
#define SIFT_CURV_THR			10		// threshold on keypoint ratio of principle curvatures,default is 10 
#define SIFT_IMG_DBL			0		// double image size before pyramid construction
#define SIFT_DESCR_WIDTH		4		// width of descriptor histogram array 
#define SIFT_DESCR_HIST_BINS	8		// number of bins per histogram in descriptor array 
#define SIFT_INIT_SIGMA			0.5		// sigmao for gaussian blur
#define SIFT_IMG_BORDER			5		// width of border in which to ignore keypoints 
#define SIFT_MAX_INTERP_STEPS	5		// maximum steps of keypoint interpolation before failure 
#define SIFT_ORI_HIST_BINS		36		// number of bins in histogram for orientation assignment 
#define SIFT_ORI_SIG_FCTR		1.5		// determines gaussian sigma for orientation assignment
#define SIFT_ORI_RADIUS			3.0 * SIFT_ORI_SIG_FCTR // the radius of the region used in orientation assignment
#define SIFT_ORI_SMOOTH_PASSES	2		// number of passes of orientation histogram smoothing 
#define SIFT_ORI_PEAK_RATIO		0.8  	// orientation magnitude relative to max that results in new feature,default 0.8
#define SIFT_DESCR_SCL_FCTR		3.0		// the size of a single descriptor orientation histogram 
#define SIFT_DESCR_MAG_THR		0.2		// threshold on magnitude of elements of descriptor vector 
#define SIFT_INT_DESCR_FCTR		512.0	// factor used to convert floating-point descriptor to unsigned char 

#define KDTREE_BBF_MAX_NN_CHKS 200
#define NN_SQ_DIST_RATIO_THR   0.7

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef MIN
#define MIN(x,y) ( ( x < y )? x : y )
#endif

#ifndef MAX
#define MAX(x,y) ( ( x > y )? x : y )
#endif
#ifndef ABS
#define ABS(x) ( ( x < 0 )? -x : x )
#endif

#ifndef INT_MAX
#define INT_MAX       2147483647    /* maximum (signed) int value */
#endif // INT_MAX

/* max feature descriptor length */
#define FEATURE_MAX_D 128


/*
	Structure to represent an affine invariant image feature.  The fields
	x, y, a, b, c represent the affine region around the feature:
	a(x-u)(x-u) + 2b(x-u)(y-v) + c(y-v)(y-v) = 1
*/


typedef struct stPoint2D64f
{
	double x;
	double y;
}
Point2D64f;


struct feature
{
	float x;                      /* x coord */
	float y;                      /* y coord */

	int d;                         /* descriptor length */
	float descr[FEATURE_MAX_D];   /* descriptor */
	
	float scl;                    /* scale of a Lowe-style feature */
	float ori;                    /* orientation of a Lowe-style feature */

	int nclass;                     /* all-purpose feature class */
	struct feature* fwd_match;     /* matching feature from forward image */
	struct feature* bck_match;     /* matching feature from backmward image */
	struct feature* mdl_match;     /* matching feature from model */
	Point2D64f img_pt;           /* location in image */
	Point2D64f mdl_pt;           /* location in model */
	void* feature_data;            /* user-definable data */
	int   index;                   //added by xdh, denote the index of feature 
};

typedef struct POS_STRUCTURE
{
	int    nImageIndex;
	double f;
	double u0,v0;
	double R[9];
	double T[3];
	double k1,k2;
	double l,r,t,b;
	double aveofdem;
	double xs,ys,zs;        //camera position
	double lon,lat;         //pos - lontitude, latitude 
	double altitude;        //pos - altitude
	double yaw,pitch,roll;  //rotation angle
	double gx,gy,gz;        //UTM coordinates
	int    zoneNumber;
}stPOS;


//parameters for absolute pose estimation
typedef struct ABS_POS_PARAMS
{
	double scale;
	double R[9];
	double T[3];
}stAbsPOS;



#endif
