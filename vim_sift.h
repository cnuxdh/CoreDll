#ifndef VIM_SIFT_H
#define VIM_SIFT_H

#include "exports.h"
#include "defs.h"


#ifdef OPENCV_1X 
	#include "cv.h"
	#include "highgui.h"
	#include "cxcore.h"
#else
	#include "opencv2/core/core.hpp"
	#include "opencv2/imgproc/imgproc_c.h"
	#include "opencv2/highgui/highgui.hpp"
	using namespace cv;
#endif




/* Interpolates a histogram peak from left, center, and right values */
#define interp_hist_peak( l, c, r ) ( 0.5 * ((l)-(r)) / ((l) - 2.0*(c) + (r)) )

/************************* Local Function Prototypes *************************/

IplImage*	ipl_create_init_img( IplImage*, int n, double dP);
IplImage*	ipl_convert2gray32( IplImage* img);
IplImage*** ipl_build_gauss_pyramids( IplImage* img, int n, int m, double dP);
IplImage*	ipl_down_sample( IplImage* img);
IplImage*** ipl_build_dog_pyramids( IplImage*** img, int n, int s);
CvSeq*		ipl_scale_space_extrema( IplImage*** , int, int, double, int, CvMemStorage*);
CvMat*		ipl_deriv_3D( IplImage***, int, int, int, int );
CvMat*		ipl_hessian_3D( IplImage***, int, int, int, int );

int is_extremum( IplImage***, int, int, int, int );
struct feature* interp_extremum( IplImage***, int, int, int, int, int, double);
void interp_step( IplImage***, int, int, int, int, double*, double*, double* );
double interp_contr( IplImage***, int, int, int, int, double, double, double );
struct feature* new_feature( void );
int is_too_edge_like( IplImage*, int, int, int );
void calc_feature_scales( CvSeq*, double, int );
void adjust_for_img_dbl( CvSeq* );
void calc_feature_oris( CvSeq*, IplImage*** );
double* ori_hist( IplImage*, int, int, int, int, double );
int calc_grad_mag_ori( IplImage*, int, int, double*, double* );
void smooth_ori_hist( double*, int );
double dominant_ori( double*, int );
void add_good_ori_features( CvSeq*, double*, int, double, struct feature* );
struct feature* clone_feature( struct feature* );
void compute_descriptors( CvSeq*, IplImage***, int, int );
double*** descr_hist( IplImage*, int, int, double, double, int, int );
void interp_hist_entry( double***, double, double, double, double, int, int);
void hist_to_descr_descriptor_vector( double***, int, int, struct feature* );

int  feature_cmpare( void*, void*);
void release_descr_hist( double****, int );
void release_pyr( IplImage****, int, int );

int DLL_EXPORT sift_features( IplImage* img, struct feature** feat );
int _sift_features( IplImage* img, struct feature** feat, int intvls,
						  double sigma, double contr_thr, int curv_thr,
						  int img_dbl, int descr_width, int descr_hist_bins );
/*
	Displays a set of features on an image
	img  - image on which to display features
	feat - array of Oxford-type features
	n    - number of features
*/
DLL_EXPORT void draw_features( IplImage* img, struct feature* feat, int n );

void draw_lowe_features( IplImage*, struct feature*, int );
void draw_lowe_feature( IplImage*, struct feature*, CvScalar );

DLL_EXPORT IplImage* extractPatch(struct feature ptFeat, IplImage* pImage, double relativeAngle=0);


/**
	A function to get a pixel value from a 32-bit floating-point image.

	img - an image
	r -  row
	c -  column
	return:
	Returns the value of the pixel at (\a r, \a c) in \a img
*/

static __inline float pixval32f(IplImage* img, int r, int c)
{
	return ((float*)(img->imageData + img->widthStep*r))[c];
}

/*-------------------------------------------------------------
Function: save sift features of one image into the file, wrote by Xie Donghai, 2012.6.27
input: 
 pFeat:    sift feature vectors, scale, orientation, etc.
 nFeat:    the number of feature points
 filepath: file to save all feature points
-------------------------------------------------------------*/
int DLL_EXPORT writeSift(struct feature* pFeat, int nFeat, char* filepath);

int DLL_EXPORT writeSiftForBundler(struct feature* pFeat, int nFeat, char* filepath);

int DLL_EXPORT writeSiftForBundlerBin(struct feature* pFeat, int nFeat, char* filepath);


/************************* END Local Functions *************************/
#endif