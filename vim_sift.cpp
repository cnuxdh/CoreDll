
/*
  vim_sift.c
  This file contains sift functions.
  
  reference: 
  1.Lowe,D.G. Distinctive image features from scale-invariant keypoints. 
  International Journal of Computer Vision,2004,2,vol.60,91-110.
  2.
*/

#include <string.h>
#include <math.h>
#include <stdio.h>

//#include "CommonFuncs.h"

#include "vim_sift.h"
#include "vim_imgs.h"


/*********************** Functions prototyped in sift.h **********************/

int sift_features( IplImage* img, struct feature** feat )
{
	return _sift_features( img, feat, SIFT_INTVLS, SIFT_SIGMA, SIFT_CONTR_THR,
							SIFT_CURV_THR, SIFT_IMG_DBL, SIFT_DESCR_WIDTH,
							SIFT_DESCR_HIST_BINS );
}


/**
	Finds SIFT features in an image using user-specified parameter values.  
	All	detected features are stored in the array pointed to by \a feat.
	Returns the number of keypoints stored in feat or -1 on failure 
*/
int _sift_features( IplImage* img, struct feature** feat, int intvls,
				   double sigma, double contr_thr, int curv_thr,
				   int img_dbl, int descr_width, int descr_hist_bins )
{
	IplImage* init_img;
	IplImage* pImage;
	IplImage*** gauss_pyramids, *** dog_pyramids;
	CvMemStorage* storage;
	CvSeq* features;
	int octvs, i, n = 0;

	/* check arguments */
	if( (! img )||( ! feat ))
		fatal_error( "pointer error");
		
	/* build scale space pyramid; smallest dimension of top level is ~4 pixels */
	init_img = ipl_create_init_img( img, img_dbl, sigma );	    
	
	//为了方便保存，将32bit单精度浮点图像转换为8bit的灰度图像,added by xdh,2009.9.17
	//pImage = cvCreateImage( cvGetSize(init_img), IPL_DEPTH_8U, 1);
	//cvConvertScale(init_img, pImage, 255, 0);
	//cvCvtColor(init_img, pImage, CV_GRAY2RGB);
	//cvSaveImage("d:\\init.bmp", pImage);

	octvs = (int) (log( (double)MIN( init_img->width, init_img->height ) ) / log(2.0) - 6);
	
	gauss_pyramids = ipl_build_gauss_pyramids( init_img, octvs, intvls, sigma );
	
	dog_pyramids = ipl_build_dog_pyramids( gauss_pyramids, octvs, intvls );

	storage = cvCreateMemStorage( 0 );
	features = ipl_scale_space_extrema( dog_pyramids, octvs, intvls, contr_thr, \
		curv_thr, storage );

	calc_feature_scales( features, sigma, intvls );
	
	if( img_dbl )
		adjust_for_img_dbl( features );
	
	calc_feature_oris( features, gauss_pyramids );
	
	compute_descriptors( features, gauss_pyramids, descr_width, descr_hist_bins );

	/* sort features by decreasing scale and move from CvSeq to array */
	cvSeqSort( features, (CvCmpFunc)feature_cmpare, NULL );
	n = features->total;
	*feat = (struct feature*)calloc( n, sizeof(struct feature) );
	*feat = (struct feature*)cvCvtSeqToArray( features, *feat, CV_WHOLE_SEQ );
	for( i = 0; i < n; i++ )
	{
		free( (*feat)[i].feature_data );
		(*feat)[i].feature_data = NULL;
	}

	cvReleaseMemStorage( &storage );
	cvReleaseImage( &init_img );
	release_pyr( &gauss_pyramids, octvs, intvls + 3 );
	release_pyr( &dog_pyramids, octvs, intvls + 2 );
	return n;
}

/* 
	Extract feature patch from the image
    
*/
/*
IplImage* extractPatch(struct feature ptFeat, IplImage* pImage, double relativeAngle)
{
	int i,j;
	int ht = pImage->height;
	int wd = pImage->width;

	//patch rectangle
	double radius = ptFeat.scl*5;
    int l = max(0, ptFeat.x - radius);
	int r = min(wd-1, ptFeat.x+radius);
	int t = max(0, ptFeat.y - radius);
	int b = min(ht-1, ptFeat.y+radius);

	IplImage* pPatch = cvCreateImage( cvSize(r-l, b-t), 8, 1 );

	for(j=t; j<b; j++)
		for(i=l; i<r; i++)
		{
			int tx = i-(l+r)*0.5;
			int ty = j-(t+b)*0.5;
			double x = tx*cos(relativeAngle) - ty*sin(relativeAngle) + (l+r)*0.5;
			double y = tx*sin(relativeAngle) + ty*cos(relativeAngle) + (t+b)*0.5;
			pPatch->imageData[ (j-t)*pPatch->widthStep+(i-l) ] = BilinearPixel(x,y,pImage);//pImage->imageData[y*pImage->widthStep+x];
		}

	return pPatch;
}
*/



/************************ Functions prototyped here **************************/
/*
	Converts an image to 8-bit grayscale and Gaussian-smooths it.

	img - input image
	img_dbl  - if true, image is doubled in size prior to smoothing
	sigma - total std of Gaussian smoothing
*/
IplImage* ipl_create_init_img( IplImage* img, int img_dbl, double sigma )
{
	IplImage* gray, * dbl;
	float sig_diff;

	//将图像从8bit灰度图像变为32bit浮点图像
	gray = ipl_convert2gray32( img );
	if( img_dbl )
	{
		sig_diff = (float) sqrt( sigma * sigma - SIFT_INIT_SIGMA * SIFT_INIT_SIGMA * 4 );
		
		dbl = cvCreateImage( cvSize( img->width*2, img->height*2 ),
			IPL_DEPTH_32F, 1 );

		//将图像进行缩放
		cvResize( gray, dbl, CV_INTER_CUBIC );
		
		//对图像进行平滑
		cvSmooth( dbl, dbl, CV_GAUSSIAN, 0, 0, sig_diff, sig_diff );
		
		cvReleaseImage( &gray );
		return dbl;
	}
	else
	{
		sig_diff = (float) sqrt( sigma * sigma - SIFT_INIT_SIGMA * SIFT_INIT_SIGMA );
		cvSmooth( gray, gray, CV_GAUSSIAN, 0, 0, sig_diff, sig_diff );
		return gray;
	}
}

/*
	Converts an image to 32-bit grayscale
	img  - a 3-channel 8-bit color (BGR) or 8-bit gray image
	return:
	Returns a 32-bit grayscale image
	added by xdh: 这意味着sift是在灰度图像上来完成的。
*/
IplImage* ipl_convert2gray32( IplImage* img )
{
	IplImage* gray8, * gray32;
	/* int r, c; */

	gray8 = cvCreateImage( cvGetSize(img), IPL_DEPTH_8U, 1 );
	gray32 = cvCreateImage( cvGetSize(img), IPL_DEPTH_32F, 1 );

	if( img->nChannels == 1 )
		gray8 = (IplImage*)cvClone( img );
	else
		cvCvtColor( img, gray8, CV_RGB2GRAY ); //将三波段图像转换为灰度图像

	//将整数转换为浮点类型
	cvConvertScale( gray8, gray32, 1.0 / 255.0, 0 );

	cvReleaseImage( &gray8 );
	return gray32;
}

/*
	Builds Gaussian scale space pyramid from an image
	base  - base image of the pyramid
	octvs  - number of octaves of scale space
	intvls  - number of intervals per octave
	sigma  - amount of Gaussian smoothing per octave
	return:
    Returns a Gaussian scale space pyramid as an octvs x (intvls + 3) array
*/
IplImage*** ipl_build_gauss_pyramids( IplImage* base, int octvs,
							int intvls, double sigma )
{
	IplImage*** gauss_pyramids;
	double* sig = (double*)calloc( intvls + 3, sizeof(double));
	double sig_total, sig_prev, k;
	int i, j;

	gauss_pyramids = (IplImage***)calloc( octvs, sizeof( IplImage** ) );
	for( i = 0; i < octvs; i++ )
		gauss_pyramids[i] = (IplImage**)calloc( intvls + 3, sizeof( IplImage* ) );
	/*
		precompute Gaussian sigmas using the following formula:
		\sigma_{total}^2 = \sigma_{i}^2 + \sigma_{i-1}^2
	*/
	sig[0] = sigma;
	//printf(" sigma: \n %lf \n", sig[0]);
	k = pow( 2.0, 1.0 / intvls );
	for( i = 1; i < intvls + 3; i++ )
	{
		sig_prev = pow( k, i - 1 ) * sigma;
		sig_total = sig_prev * k;
		sig[i] = sqrt( sig_total * sig_total - sig_prev * sig_prev );
		//printf("%lf \n", sig[i]);
	}    

	for( j = 0; j < octvs; j++ )
		for( i = 0; i < intvls + 3; i++ )
		{
			if( j == 0  &&  i == 0 )
				gauss_pyramids[j][i] = cvCloneImage(base);

			/* base of new octvave is halved image from end of previous octave */
			else if( i == 0 )
				gauss_pyramids[j][i] = ipl_down_sample( gauss_pyramids[j-1][intvls] );

			/* blur the current octave's last image to create the next one */
			else
			{
				gauss_pyramids[j][i] = cvCreateImage( cvGetSize(gauss_pyramids[j][i-1]),
					IPL_DEPTH_32F, 1 );
				cvSmooth( gauss_pyramids[j][i-1], gauss_pyramids[j][i],
					CV_GAUSSIAN, 0, 0, sig[i], sig[i] );
			}
		}
	free( sig );
	return gauss_pyramids;
}

/*
	Downsamples an image to a quarter of its size (half in each dimension)
	using nearest-neighbor interpolation
	img an image
	return Returns an image whose dimensions are half those of img
*/
IplImage* ipl_down_sample( IplImage* img )
{
	IplImage* smaller = cvCreateImage( cvSize(img->width / 2, img->height / 2),
		img->depth, img->nChannels );
	cvResize( img, smaller, CV_INTER_NN );

	return smaller;
}

/*
	Builds a difference of Gaussians scale space pyramid by subtracting adjacent
	intervals of a Gaussian pyramid

	gauss_pyramids - Gaussian scale-space pyramid
	octvs  - number of octaves of scale space
	intvls - number of intervals per octave
+-
	return:
	Returns a difference of Gaussians scale space pyramid as an octvs x (intvls + 2) array
*/
IplImage*** ipl_build_dog_pyramids( IplImage*** gauss_pyramids, int octvs, int intvls )
{
	IplImage*** dog_pyramids;
	int i, j;

	dog_pyramids = (IplImage***)calloc( octvs, sizeof( IplImage** ) );
	for( i = 0; i < octvs; i++ )
		dog_pyramids[i] = (IplImage**)calloc( intvls + 2, sizeof(IplImage*) );

	for( j = 0; j < octvs; j++ )
		for( i = 0; i < intvls + 2; i++ )
		{
			dog_pyramids[j][i] = cvCreateImage( cvGetSize(gauss_pyramids[j][i]),
				IPL_DEPTH_32F, 1 );
			cvSub( gauss_pyramids[j][i+1], gauss_pyramids[j][i], dog_pyramids[j][i], NULL );
		}
		

	return dog_pyramids;
}

/*
	Determines whether a pixel is a scale-space extremum by comparing it to it's
	3x3x3 pixel neighborhood.

	dog_pyramids -  DoG scale space pyramid
	octv  - pixel's scale space octave
	intvl -  pixel's within-octave interval
	r -  pixel's image row
	c - pixel's image col

	return:
	Returns 1 if the specified pixel is an extremum (max or min) among
		it's 3x3x3 pixel neighborhood.
*/
int is_extremum( IplImage*** dog_pyramids, int octv, int intvl, int r, int c )
{
	float val = pixval32f( dog_pyramids[octv][intvl], r, c );
	int i, j, k;

	/* check for maximum */
	if( val > 0 )
	{
		for( i = -1; i <= 1; i++ )
		{	
			for( j = -1; j <= 1; j++ )
			{
				for( k = -1; k <= 1; k++ )
				{
					if( val < pixval32f( dog_pyramids[octv][intvl+i], r + j, c + k ) )
						return 0;
				}
			}
		}
	}

	/* check for minimum */
	else
	{
		for( i = -1; i <= 1; i++ )
			for( j = -1; j <= 1; j++ )
				for( k = -1; k <= 1; k++ )
					if( val > pixval32f( dog_pyramids[octv][intvl+i], r + j, c + k ) )
						return 0;
	}

	return 1;
}

/*
	Interpolates a scale-space extremum's location and scale to subpixel 
	accuracy to form an image feature. Rejects features with low contrast.
	Based on Section 4 of Lowe's paper.  

	dog_pyramids -  DoG scale space pyramid
	octv -  feature's octave of scale space
	intvl -  feature's within-octave interval
	r -  feature's image row
	c -  feature's image column
	intvls  - total intervals per octave
	contr_thr -  threshold on feature contrast

	return:
    Returns the feature resulting from interpolation of the given
	    parameters or NULL if the given location could not be interpolated or
		if contrast at the interpolated loation was too low.  
	If a feature is returned, its scale, orientation, and descriptor are 
	    yet to be determined.
*/
struct feature* interp_extremum( IplImage*** dog_pyramids, int octv, int intvl,
								int r, int c, int intvls, double contr_thr )
{
	struct feature* feat;
	struct detection_data* ddata;
	double xi, xr, xc, contr;
	int i = 0;

	while( i < SIFT_MAX_INTERP_STEPS )
	{
		interp_step( dog_pyramids, octv, intvl, r, c, &xi, &xr, &xc );
		if( ABS( xi ) < 0.5  &&  ABS( xr ) < 0.5  &&  ABS( xc ) < 0.5 )
			break;

		c += cvRound( xc );
		r += cvRound( xr );
		intvl += cvRound( xi );

		if( intvl < 1  || intvl > intvls  ||c < SIFT_IMG_BORDER  || \
			r < SIFT_IMG_BORDER  ||
			c >= dog_pyramids[octv][0]->width - SIFT_IMG_BORDER  ||
			r >= dog_pyramids[octv][0]->height - SIFT_IMG_BORDER )
		{
			return NULL;
		}

		i++;
	}

	/* ensure convergence of interpolation */
	if( i >= SIFT_MAX_INTERP_STEPS )
		return NULL;

	contr = interp_contr( dog_pyramids, octv, intvl, r, c, xi, xr, xc );
	if( ABS( contr ) < contr_thr / intvls )
		return NULL;

	feat = new_feature();
	ddata = feat_detection_data( feat );
	feat->img_pt.x = feat->x = ( c + xc ) * pow( 2.0, octv );
	feat->img_pt.y = feat->y = ( r + xr ) * pow( 2.0, octv );
	ddata->r = r;
	ddata->c = c;
	ddata->octv = octv;
	ddata->intvl = intvl;
	ddata->subintvl = xi;

	return feat;
}

/*
	Performs one step of extremum interpolation.  
	Based on Eqn. (3) in Lowe's paper.

	dog_pyramids  - difference of Gaussians scale space pyramid
	octv  - octave of scale space
	intvl  - interval being interpolated
	r -  row being interpolated
	c -  column being interpolated
	xi -  output as interpolated subpixel increment to interval
	xr -  output as interpolated subpixel increment to row
	xc -  output as interpolated subpixel increment to col
*/

void interp_step( IplImage*** dog_pyramids, int octv, int intvl, int r, int c,
				 double* xi, double* xr, double* xc )
{
	CvMat* dD, * H, * H_inv, X;
	double x[3] = { 0 };

	dD = ipl_deriv_3D( dog_pyramids, octv, intvl, r, c );
	H = ipl_hessian_3D( dog_pyramids, octv, intvl, r, c );
	H_inv = cvCreateMat( 3, 3, CV_64FC1 );
	cvInvert( H, H_inv, CV_SVD );
	cvInitMatHeader( &X, 3, 1, CV_64FC1, x, CV_AUTOSTEP );
	cvGEMM( H_inv, dD, -1, NULL, 0, &X, 0 );

	cvReleaseMat( &dD );
	cvReleaseMat( &H );
	cvReleaseMat( &H_inv );

	*xi = x[2];
	*xr = x[1];
	*xc = x[0];
}

/*
	Computes the partial derivatives in x, y, and scale of a pixel in the DoG
	scale space pyramid.

	dog_pyramids -  DoG scale space pyramid
	octv -  pixel's octave in dog_pyramids
	intvl -  pixel's interval in octv
	r -  pixel's image row
	c - pixel's image col

	return:
	Returns the vector of partial derivatives for pixel I
		{ dI/dx, dI/dy, dI/ds }^T as a CvMat*
*/
CvMat* ipl_deriv_3D( IplImage*** dog_pyramids, int octv, int intvl, int r, int c )
{
	CvMat* dI;
	double dx, dy, ds;

	dx = ( pixval32f( dog_pyramids[octv][intvl], r, c+1 ) -
		pixval32f( dog_pyramids[octv][intvl], r, c-1 ) ) / 2.0;
	dy = ( pixval32f( dog_pyramids[octv][intvl], r+1, c ) -
		pixval32f( dog_pyramids[octv][intvl], r-1, c ) ) / 2.0;
	ds = ( pixval32f( dog_pyramids[octv][intvl+1], r, c ) -
		pixval32f( dog_pyramids[octv][intvl-1], r, c ) ) / 2.0;

	dI = cvCreateMat( 3, 1, CV_64FC1 );
	cvmSet( dI, 0, 0, dx );
	cvmSet( dI, 1, 0, dy );
	cvmSet( dI, 2, 0, ds );

	return dI;
}

/*
	Computes the 3D Hessian matrix for a pixel in the DoG scale space pyramid.

	dog_pyramids -  DoG scale space pyramid
	octv -  pixel's octave in dog_pyramids
	intvl -  pixel's interval in octv
	r -  pixel's image row
	c -  pixel's image col

	return:
	Returns the Hessian matrix (below) for pixel I as a CvMat*

		/ Ixx  Ixy  Ixs \ <BR>
		| Ixy  Iyy  Iys | <BR>
		\ Ixs  Iys  Iss /
*/
CvMat* ipl_hessian_3D( IplImage*** dog_pyramids, int octv, int intvl, int r, int c )
{
	CvMat* H;
	double v, dxx, dyy, dss, dxy, dxs, dys;

	v = pixval32f( dog_pyramids[octv][intvl], r, c );
	dxx = ( pixval32f( dog_pyramids[octv][intvl], r, c+1 ) + 
			pixval32f( dog_pyramids[octv][intvl], r, c-1 ) - 2 * v );
	dyy = ( pixval32f( dog_pyramids[octv][intvl], r+1, c ) +
			pixval32f( dog_pyramids[octv][intvl], r-1, c ) - 2 * v );
	dss = ( pixval32f( dog_pyramids[octv][intvl+1], r, c ) +
			pixval32f( dog_pyramids[octv][intvl-1], r, c ) - 2 * v );
	dxy = ( pixval32f( dog_pyramids[octv][intvl], r+1, c+1 ) -
			pixval32f( dog_pyramids[octv][intvl], r+1, c-1 ) -
			pixval32f( dog_pyramids[octv][intvl], r-1, c+1 ) +
			pixval32f( dog_pyramids[octv][intvl], r-1, c-1 ) ) / 4.0;
	dxs = ( pixval32f( dog_pyramids[octv][intvl+1], r, c+1 ) -
			pixval32f( dog_pyramids[octv][intvl+1], r, c-1 ) -
			pixval32f( dog_pyramids[octv][intvl-1], r, c+1 ) +
			pixval32f( dog_pyramids[octv][intvl-1], r, c-1 ) ) / 4.0;
	dys = ( pixval32f( dog_pyramids[octv][intvl+1], r+1, c ) -
			pixval32f( dog_pyramids[octv][intvl+1], r-1, c ) -
			pixval32f( dog_pyramids[octv][intvl-1], r+1, c ) +
			pixval32f( dog_pyramids[octv][intvl-1], r-1, c ) ) / 4.0;

	H = cvCreateMat( 3, 3, CV_64FC1 );
	cvmSet( H, 0, 0, dxx );
	cvmSet( H, 0, 1, dxy );
	cvmSet( H, 0, 2, dxs );
	cvmSet( H, 1, 0, dxy );
	cvmSet( H, 1, 1, dyy );
	cvmSet( H, 1, 2, dys );
	cvmSet( H, 2, 0, dxs );
	cvmSet( H, 2, 1, dys );
	cvmSet( H, 2, 2, dss );

	return H;
}

/*
	Calculates interpolated pixel contrast.  
	Based on Eqn. (3) in Lowe's paper.

	dog_pyramids  - difference of Gaussians scale space pyramid
	octv  - octave of scale space
	intvl -  within-octave interval
	r -  pixel row
	c  - pixel column
	xi -  interpolated subpixel increment to interval
	xr -  interpolated subpixel increment to row
	xc -  interpolated subpixel increment to col

	return:
	 Returns interpolated contrast.
*/
double interp_contr( IplImage*** dog_pyramids, int octv, int intvl, int r,
					int c, double xi, double xr, double xc )
{
	CvMat* dD, X, T;
	double t[1], x[3] = { xc, xr, xi };

	cvInitMatHeader( &X, 3, 1, CV_64FC1, x, CV_AUTOSTEP );
	cvInitMatHeader( &T, 1, 1, CV_64FC1, t, CV_AUTOSTEP );
	dD = ipl_deriv_3D( dog_pyramids, octv, intvl, r, c );
	cvGEMM( dD, &X, 1, NULL, 0, &T,  CV_GEMM_A_T );
	cvReleaseMat( &dD );

	return pixval32f( dog_pyramids[octv][intvl], r, c ) + t[0] * 0.5;
}

/*
	Determines whether a feature is too edge like to be stable by computing the
	ratio of principal curvatures at that feature.  
	Based on Section 4.1 of Lowe's paper.

	dog_img  - image from the DoG pyramid in which feature was detected
	r  - feature row
	c  - feature col
	curv_thr  - high threshold on ratio of principal curvatures

	return:
	Returns 0 if the feature at (r,c) in dog_img is sufficiently corner-like or 1 otherwise.
*/
int is_too_edge_like( IplImage* dog_img, int r, int c, int curv_thr )
{
	double d, dxx, dyy, dxy, tr, det;

	/* principal curvatures are computed using the trace and det of Hessian */
	d = pixval32f(dog_img, r, c);
	dxx = pixval32f( dog_img, r, c+1 ) + pixval32f( dog_img, r, c-1 ) - 2 * d;
	dyy = pixval32f( dog_img, r+1, c ) + pixval32f( dog_img, r-1, c ) - 2 * d;
	dxy = ( pixval32f(dog_img, r+1, c+1) - pixval32f(dog_img, r+1, c-1) -
			pixval32f(dog_img, r-1, c+1) + pixval32f(dog_img, r-1, c-1) ) / 4.0;
	tr = dxx + dyy;
	det = dxx * dyy - dxy * dxy;

	/* negative determinant -> curvatures have different signs; reject feature */
	if( det <= 0 )
		return 1;

	if( tr * tr / det < ( curv_thr + 1.0 )*( curv_thr + 1.0 ) / curv_thr )
		return 0;

	return 1;
}

/*
	Calculates characteristic scale for each feature in an array.

	features  - array of features
	sigma  - amount of Gaussian smoothing per octave of scale space
	intvls  - intervals per octave of scale space
*/
void calc_feature_scales( CvSeq* features, double sigma, int intvls )
{
	struct feature* feat;
	struct detection_data* ddata;
	double intvl;
	int i, n;

	n = features->total;
	for( i = 0; i < n; i++ )
	{
		feat = CV_GET_SEQ_ELEM( struct feature, features, i );
		ddata = feat_detection_data( feat );
		intvl = ddata->intvl + ddata->subintvl;
		feat->scl = sigma * pow( 2.0, ddata->octv + intvl / intvls );
		ddata->scl_octv = sigma * pow( 2.0, intvl / intvls );
	}
}

/*
	Halves feature coordinates and scale in case the input image was doubled
	prior to scale space construction.

	features  - array of features
*/
void adjust_for_img_dbl( CvSeq* features )
{
	struct feature* feat;
	int i, n;

	n = features->total;
	for( i = 0; i < n; i++ )
	{
		feat = CV_GET_SEQ_ELEM( struct feature, features, i );
		feat->x /= 2.0;
		feat->y /= 2.0;
		feat->scl /= 2.0;
		feat->img_pt.x /= 2.0;
		feat->img_pt.y /= 2.0;
	}
}

/*
	Computes a canonical orientation for each image feature in an array.  
	Based on Section 5 of Lowe's paper.  
	This function adds features to the array when there is more than 
	one dominant orientation at a given feature location.

	features  - an array of image features
	gauss_pyramids  - Gaussian scale space pyramid
*/
void calc_feature_oris( CvSeq* features, IplImage*** gauss_pyramids )
{
	struct feature* feat;
	struct detection_data* ddata;
	double* hist;
	double omax;
	int i, j, n = features->total;

	for( i = 0; i < n; i++ )
	{
		feat = (struct feature*)malloc( sizeof( struct feature ) );
		cvSeqPopFront( features, feat );
		ddata = feat_detection_data( feat );
		hist = ori_hist( gauss_pyramids[ddata->octv][ddata->intvl],
						ddata->r, ddata->c, SIFT_ORI_HIST_BINS,
						cvRound( SIFT_ORI_RADIUS * ddata->scl_octv ),
						SIFT_ORI_SIG_FCTR * ddata->scl_octv );
		for( j = 0; j < SIFT_ORI_SMOOTH_PASSES; j++ )
			smooth_ori_hist( hist, SIFT_ORI_HIST_BINS );
		omax = dominant_ori( hist, SIFT_ORI_HIST_BINS );
		add_good_ori_features( features, hist, SIFT_ORI_HIST_BINS,
								omax * SIFT_ORI_PEAK_RATIO, feat );
		free( ddata );
		free( feat );
		free( hist );
	}
}

/*
	Computes a gradient orientation histogram at a specified pixel.

	img  - image
	r  -  pixel row
	c  - pixel col
	n  - number of histogram bins
	rad   - radius of region over which histogram is computed
	sigma -  std for Gaussian weighting of histogram entries

	return:
	 Returns an n-element array containing an orientation histogram
		representing orientations between 0 and 2 PI.
*/
double* ori_hist( IplImage* img, int r, int c, int n, int rad, double sigma)
{
	double* hist;
	double mag, ori, w, exp_denom, PI2 = CV_PI * 2.0;
	int bin, i, j;

	hist = (double*)calloc( n, sizeof( double ) );
	exp_denom = 2.0 * sigma * sigma;
	for( i = -rad; i <= rad; i++ )
		for( j = -rad; j <= rad; j++ )
			if( calc_grad_mag_ori( img, r + i, c + j, &mag, &ori ) )
			{
				//printf("mag: %lf ori: %lf \n", mag, ori);
				w = exp( -( i*i + j*j ) / exp_denom );
				bin = cvRound( n * ( ori + CV_PI ) / PI2 );
				bin = ( bin < n )? bin : 0;
				hist[bin] += w * mag;
			}

	return hist;
}

/*
	Calculates the gradient magnitude and orientation at a given pixel.

	img - image
	r - pixel row
	c - pixel col
	mag - output as gradient magnitude at pixel (r,c)
	ori - output as gradient orientation at pixel (r,c)

	return:
	 Returns 1 if the specified pixel is a valid one and sets mag and
		ori accordingly; otherwise returns 0
*/
int calc_grad_mag_ori( IplImage* img, int r, int c, double* mag, double* ori )
{
	double dx, dy;

	if( r > 0  &&  r < img->height - 1  &&  c > 0  &&  c < img->width - 1 )
	{
		dx = pixval32f( img, r, c+1 ) - pixval32f( img, r, c-1 );
		dy = pixval32f( img, r-1, c ) - pixval32f( img, r+1, c ); //y direction from bottom to top, xdh, 2013.4.30
		*mag = sqrt( dx*dx + dy*dy );
		*ori = atan2( dy, dx );
		return 1;
	}

	else
		return 0;
}

/*
	Adds features to an array for every orientation in a histogram greater than
	a specified threshold.

	features  - new features are added to the end of this array
	hist  - orientation histogram
	n  - number of bins in hist
	mag_thr  - new features are added for entries in hist greater than this
	feat new features are clones of this with different orientations
*/
void add_good_ori_features( CvSeq* features, double* hist, int n,
						   double mag_thr, struct feature* feat )
{
	struct feature* new_feat;
	double bin, PI2 = CV_PI * 2.0;
	int l, r, i;

	for( i = 0; i < n; i++ )
	{
		l = ( i == 0 )? n - 1 : i-1;
		r = ( i + 1 ) % n;

		if( hist[i] > hist[l]  &&  hist[i] > hist[r]  &&  hist[i] >= mag_thr )
		{
			bin = i + interp_hist_peak( hist[l], hist[i], hist[r] );
			bin = ( bin < 0 )? n + bin : ( bin >= n )? bin - n : bin;
			new_feat = clone_feature( feat );
			new_feat->ori = ( ( PI2 * bin ) / n ) - CV_PI;
			cvSeqPush( features, new_feat );
			free( new_feat );
		}
	}
}

/*
	Computes feature descriptors for features in an array.  
	Based on Section 6 of Lowe's paper.

	features - array of features
	gauss_pyramids - Gaussian scale space pyramid
	d - width of 2D array of orientation histograms
	n - number of bins per orientation histogram
*/
void compute_descriptors( CvSeq* features, IplImage*** gauss_pyramids, int d, int n)
{
	struct feature* feat;
	struct detection_data* ddata;
	double*** hist;
	int i, k = features->total;

	for( i = 0; i < k; i++ )
	{
		feat = CV_GET_SEQ_ELEM( struct feature, features, i );
		ddata = feat_detection_data( feat );
		hist = descr_hist( gauss_pyramids[ddata->octv][ddata->intvl], ddata->r,
			ddata->c, feat->ori, ddata->scl_octv, d, n );
		hist_to_descr_descriptor_vector( hist, d, n, feat );
		release_descr_hist( &hist, d );
	}
}

/*
	Computes the 2D array of orientation histograms that form the feature
	descriptor.  Based on Section 6.1 of Lowe's paper.

	img - image used in descriptor computation
	r - row coord of center of orientation histogram array
	c - column coord of center of orientation histogram array
	ori - canonical orientation of feature whose descr is being computed
	scl  - scale relative to img of feature whose descr is being computed
	d - width of 2d array of orientation histograms
	n - bins per orientation histogram
ó
	return:
	Returns a d x d array of n-bin orientation histograms.
*/
double*** descr_hist( IplImage* img, int r, int c, double ori, double scl, int d, int n )
{
	double*** hist;
	double cos_t, sin_t, hist_width, exp_denom, r_rot, c_rot, grad_mag,
		grad_ori, w, rbin, cbin, obin, bins_per_rad, PI2 = 2.0 * CV_PI;
	int radius, i, j;

	hist = (double***)calloc( d, sizeof( double** ) );
	for( i = 0; i < d; i++ )
	{
		hist[i] = (double**)calloc( d, sizeof( double* ) );
		for( j = 0; j < d; j++ )
			hist[i][j] = (double*)calloc( n, sizeof( double ) );
	}

	cos_t = cos( ori );
	sin_t = sin( ori );
	bins_per_rad = n / PI2;
	exp_denom = d * d * 0.5;
	hist_width = SIFT_DESCR_SCL_FCTR * scl;
	radius = (int) (hist_width * sqrt(2.0) * ( d + 1.0 ) * 0.5 + 0.5);
	for( i = -radius; i <= radius; i++ )
		for( j = -radius; j <= radius; j++ )
		{
			/*
			Calculate sample's histogram array coords rotated relative to ori.
			Subtract 0.5 so samples that fall e.g. in the center of row 1 (i.e.
			r_rot = 1.5) have full weight placed in row 1 after interpolation.
			*/
			c_rot = ( j * cos_t - i * sin_t ) / hist_width;
			r_rot = ( j * sin_t + i * cos_t ) / hist_width;
			rbin = r_rot + d / 2 - 0.5;
			cbin = c_rot + d / 2 - 0.5;

			if( rbin > -1.0  &&  rbin < d  &&  cbin > -1.0  &&  cbin < d )
				if( calc_grad_mag_ori( img, r + i, c + j, &grad_mag, &grad_ori ))
				{
					grad_ori -= ori;
					while( grad_ori < 0.0 )
						grad_ori += PI2;
					while( grad_ori >= PI2 )
						grad_ori -= PI2;

					obin = grad_ori * bins_per_rad;
					w = exp( -(c_rot * c_rot + r_rot * r_rot) / exp_denom );
					interp_hist_entry( hist, rbin, cbin, obin, grad_mag * w, d, n );
				}
		}

	return hist;
}

/*
	De-allocates memory held by a scale space pyramid

	pyr - scale space pyramid
	octvs - number of octaves of scale space
	n  - number of images per octave
*/
void release_pyr( IplImage**** pyr, int octvs, int n )
{
	int i, j;
	for( i = 0; i < octvs; i++ )
	{
		for( j = 0; j < n; j++ )
			cvReleaseImage( &(*pyr)[i][j] );
		free( (*pyr)[i] );
	}
	free( *pyr );
	*pyr = NULL;
}

/*
	Detects features at extrema in DoG scale space.  
	Bad features are discarded based on contrast and ratio of principal curvatures.

	dog_pyramids -  DoG scale space pyramid
	octvs  - octaves of scale space represented by dog_pyramids
	intvls  - intervals per octave
	contr_thr  - low threshold on feature contrast
	curv_thr  - high threshold on feature ratio of principal curvatures
	storage  - memory storage in which to store detected features

	return:
	Returns an array of detected features whose scales, orientations,and descriptors are yet to be determined.
*/
CvSeq* ipl_scale_space_extrema( IplImage*** dog_pyramids, int octvs, int intvls,
						   double contr_thr, int curv_thr,
						   CvMemStorage* storage )
{
	CvSeq* features;
	double prelim_contr_thr = 0.5 * contr_thr / intvls;
	struct feature* feat;
	struct detection_data* ddata;
	int j, i , r, c/*, w, h*/;

	features = cvCreateSeq( 0, sizeof(CvSeq), sizeof(struct feature), storage );
	for( j=0; j<octvs; j++ )
	{
		for( i=1; i<=intvls; i++ )
		{
			for(r = SIFT_IMG_BORDER; r < dog_pyramids[j][0]->height-SIFT_IMG_BORDER; r++)
			{
				for(c = SIFT_IMG_BORDER; c < dog_pyramids[j][0]->width-SIFT_IMG_BORDER; c++)
				{
					/* perform preliminary check on contrast */
					if( ABS( pixval32f( dog_pyramids[j][i], r, c ) ) > prelim_contr_thr )
					{
						if( is_extremum( dog_pyramids, j, i, r, c ) )
						{
							feat = interp_extremum(dog_pyramids, j, i, r, c, intvls, contr_thr);
							if( feat )
							{
								ddata = feat_detection_data( feat );
								if( ! is_too_edge_like( dog_pyramids[ddata->octv][ddata->intvl],ddata->r, ddata->c, curv_thr ) )
								{
									cvSeqPush( features, feat );
								}
								else
									free( ddata );
								free( feat );
							}
						}
					}
				}
			}
		}
	}
	return features;
}

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

/*
	Draws a single Lowe-type feature

	img - image on which to draw
	feat -  feature to be drawn
	color -  color in which to draw
*/
void draw_lowe_feature( IplImage* img, struct feature* feat, CvScalar color )
{
	int len, hlen, blen, start_x, start_y, end_x, end_y, h1_x, h1_y, h2_x, h2_y;
	double scl, ori;
	double scale = 5.0;
	double hscale = 0.75;
	CvPoint start, end, h1, h2;

	/* compute points for an arrow scaled and rotated by feat's scl and ori */
	start_x = cvRound( feat->x );
	start_y = cvRound( feat->y );
	scl = feat->scl;
	ori = feat->ori;
	len = cvRound( scl * scale );
	hlen = cvRound( scl * hscale );
	blen = len - hlen;
	end_x = cvRound( len *  cos( ori ) ) + start_x;
	end_y = cvRound( len * -sin( ori ) ) + start_y;
	h1_x = cvRound( blen *  cos( ori + CV_PI / 18.0 ) ) + start_x;
	h1_y = cvRound( blen * -sin( ori + CV_PI / 18.0 ) ) + start_y;
	h2_x = cvRound( blen *  cos( ori - CV_PI / 18.0 ) ) + start_x;
	h2_y = cvRound( blen * -sin( ori - CV_PI / 18.0 ) ) + start_y;
	start = cvPoint( start_x, start_y );
	end = cvPoint( end_x, end_y );
	h1 = cvPoint( h1_x, h1_y );
	h2 = cvPoint( h2_x, h2_y );

	cvLine( img, start, end, color, 1, 8, 0 );
	cvLine( img, end, h1, color, 1, 8, 0 );
	cvLine( img, end, h2, color, 1, 8, 0 );

	//CvFont font;	
	//cvInitFont( &font, CV_FONT_VECTOR0,0.5, 0.5, 0, 1, 8);	
}

/*
	Draws Lowe-type features

	img - image on which to draw features
	feat - array of Oxford-type features
	n - number of features
*/
void draw_lowe_features( IplImage* img, struct feature* feat, int n )
{
	CvScalar color = CV_RGB( 255, 255, 255 );
	int i;

	if( img-> nChannels > 1 )
		color = FEATURE_LOWE_COLOR;

	for( i = 0; i < n; i++ )
		draw_lowe_feature( img, feat + i, color );
}



/* -------------------------------------
Draws a set of features on an image
img  - image on which to draw features
feat - array of Oxford-type features
n - number of features
-------------------------------------*/
void draw_features( IplImage* img, struct feature* feat, int n )
{
	if( n <= 0  ||  ! feat )
		return;
	draw_lowe_features( img, feat, n );
}


int writeSift(struct feature* pFeat, int nFeat, char* filepath)
{
	int i,j;
	FILE* fp = NULL;
	int ndim;
	
	ndim = 128;
	
	fp = fopen(filepath, "w");
	if(fp==NULL)
		return 0;

	fprintf(fp, "%d %d \n", nFeat, ndim);
	for(i=0; i<nFeat; i++)
	{
		//x y scale orientation
		fprintf(fp, "%d   %lf %lf %lf %lf \n", i, pFeat[i].x, pFeat[i].y, pFeat[i].scl, pFeat[i].ori);
		for(j=0; j<ndim; j++)
		{
			//fprintf(fp, "%6.4lf ", pFeat[i].descr[j]);
			fprintf(fp, "%d ", (int)(pFeat[i].descr[j]) );
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	return 1;
}
/*
   write sift feature points into file based on bundler format, ASCII format
*/
int writeSiftForBundler(struct feature* pFeat, int nFeat, char* filepath)
{
	int i,j;
	FILE* fp = NULL;
	int ndim;

	ndim = 128;

	fp = fopen(filepath, "w");
	if(fp==NULL)
		return 0;

	fprintf(fp, "%d %d \n", nFeat, ndim);
	for(i=0; i<nFeat; i++)
	{
		//y x scale orientation
		fprintf(fp, "%lf %lf %lf %lf \n", pFeat[i].y, pFeat[i].x, pFeat[i].scl, pFeat[i].ori);
				
		int index = 0;
		for(int m=0; m<6; m++)
		{
			for(int n=0; n<20; n++)
			{				
				fprintf(fp, "%d ", (int)(pFeat[i].descr[index]) );
				index ++;
			}
			fprintf(fp, "\n");
		}
		for(int m=0; m<8; m++)
		{
			fprintf(fp,"%d ", (int)(pFeat[i].descr[index]));
			index++;
		}		
		fprintf(fp, "\n");
	}
	fclose(fp);

	return 1;

}

//save sift feature using Binary format
int writeSiftForBundlerBin(struct feature* pFeat, int nFeat, char* filepath)
{
	FILE* fp = fopen(filepath, "wb");
	if(fp == NULL)
		return 0;

	//number of feature points
	int numFeat = nFeat;
	fwrite(&nFeat, sizeof(int), 1, fp );	
	//the dim of feature 
	int ndim = 128;
	fwrite(&ndim, sizeof(int), 1, fp);    
	//write feature vector one by one
	for(int i=0; i<nFeat; i++)
	{
		//y x scale orientation
		float y,x,scale,ori;
		x = pFeat[i].x;
		y = pFeat[i].y;
		scale = pFeat[i].scl;
		ori = pFeat[i].ori;
		fwrite(&y, sizeof(float), 1, fp);
		fwrite(&x, sizeof(float), 1, fp);
		fwrite(&scale, sizeof(float), 1, fp);
		fwrite(&ori, sizeof(float), 1, fp);
		//feature vector
    fwrite(pFeat[i].descr, sizeof(float), 128, fp);
	}
	fclose(fp);

	return 1;
}