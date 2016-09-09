
#include "VideoIndex.h"

#include "../../Corelib/Hist.h"




/*
Converts a BGR image to HSV colorspace

@param bgr image to be converted

@return Returns bgr converted to a 3-channel, 32-bit HSV image with
S and V values in the range [0,1] and H value in the range [0,360]
*/
IplImage* bgr2hsv( IplImage* bgr )
{
	IplImage* bgr32f, * hsv;

	bgr32f = cvCreateImage( cvGetSize(bgr), IPL_DEPTH_32F, 3 );
	hsv = cvCreateImage( cvGetSize(bgr), IPL_DEPTH_32F, 3 );
	cvConvertScale( bgr, bgr32f, 1.0 / 255.0, 0 );
	cvCvtColor( bgr32f, hsv, CV_BGR2HSV );
	cvReleaseImage( &bgr32f );
	return hsv;
}


/*
Calculates the histogram bin into which an HSV entry falls

@param h Hue
@param s Saturation
@param v Value

@return Returns the bin index corresponding to the HSV color defined by
\a h, \a s, and \a v.
*/
/*
int histo_bin( float h, float s, float v )
{
	int hd, sd, vd;
	
	vd = MIN( (int)(v * NV / V_MAX), NV-1 );
	if( s < S_THRESH  ||  v < V_THRESH )
		return NH * NS + vd;

	
	hd = MIN( (int)(h * NH / H_MAX), NH-1 );
	sd = MIN( (int)(s * NS / S_MAX), NS-1 );
	return sd * NH + hd;
}
*/

/*
Calculates a cumulative histogram as defined above for a given array
of images

@param img an array of images over which to compute a cumulative histogram;
each must have been converted to HSV colorspace using bgr2hsv()
@param n the number of images in imgs

@return Returns an un-normalized HSV histogram for \a imgs
*/
void calc_histogram( IplImage* img, histogram* histo )
{
	IplImage* h, * s, * v;
	float* hist;
	int i, r, c, bin;

	histo->n = NH*NS + NV;
	hist = histo->histo;
	memset( hist, 0, histo->n * sizeof(float) );

	/* extract individual HSV planes from image */
	h = cvCreateImage( cvGetSize(img), IPL_DEPTH_32F, 1 );
	s = cvCreateImage( cvGetSize(img), IPL_DEPTH_32F, 1 );
	v = cvCreateImage( cvGetSize(img), IPL_DEPTH_32F, 1 );
	cvCvtPixToPlane( img, h, s, v, NULL );

	/* increment appropriate histogram bin for each pixel */
	for( r = 0; r < img->height; r++ )
		for( c = 0; c < img->width; c++ )
		{
			bin = histo_bin( /*pixval32f( h, r, c )*/((float*)(h->imageData + h->widthStep*r) )[c],
				((float*)(s->imageData + s->widthStep*r) )[c],
				((float*)(v->imageData + v->widthStep*r) )[c] );
			hist[bin] += 1;
		}
	cvReleaseImage( &h );
	cvReleaseImage( &s );
	cvReleaseImage( &v );
}


/* calculate the gray image histogram

*/
void CalculateGrayHistogram(IplImage* gray, histogram* pHist)
{	
	int ht,wd;
	int r,c;
	int nstep = 256 / pHist->n;

	memset( pHist->histo, 0, sizeof(double)*pHist->n );

	for( r = 0; r < gray->height; r++ )
		for( c = 0; c < gray->width; c++ )
		{
			int bin = (unsigned char)(gray->imageData[r*gray->widthStep+c]) / nstep ;
			pHist->histo[bin] += 1;
		}
}


/*
Normalizes a histogram so all bins sum to 1.0

@param histo a histogram
*/
void normalize_histogram( histogram* histo )
{
	float* hist;
	float sum = 0, inv_sum;
	int i, n;

	hist = histo->histo;
	n = histo->n;

	/* compute sum of all bins and multiply each bin by the sum's inverse */
	for( i = 0; i < n; i++ )
		sum += hist[i];
	inv_sum = 1.0 / sum;
	for( i = 0; i < n; i++ )
		hist[i] *= inv_sum;
}

/*
Computes squared distance metric based on the Battacharyya similarity
coefficient between histograms.

@param h1 first histogram; should be normalized
@param h2 second histogram; should be normalized

@return Returns a squared distance based on the Battacharyya similarity
coefficient between \a h1 and \a h2
*/
float histo_dist_sq( histogram* h1, histogram* h2 )
{
	float* hist1, * hist2;
	float sum = 0;
	int i, n;

	n = h1->n;
	hist1 = h1->histo;
	hist2 = h2->histo;

	/*
	According the the Battacharyya similarity coefficient,

	D = \sqrt{ 1 - \sum_1^n{ \sqrt{ h_1(i) * h_2(i) } } }
	*/
	for( i = 0; i < n; i++ )
		sum += sqrt( hist1[i]*hist2[i] );
	return 1.0 - sum;
}


int VideoSearch(IplImage* pFrame, histogram* pDstHist, CvRect* pRect)
{     
	//load video frame image and convert into HSV image
	double ratio = 1; //(double)(pFrame->height) / (double)(pFrame->width);
	int ht = pFrame->height;
	int wd = pFrame->width;
	double zs = 1.3;
	int    objWd = 32;
	int    objHt = objWd*ratio;
	int    level = log(wd/(objWd*4.0)) / log(zs);//log2( (wd/32) );
	int    i,j,k;
	int    nwindow = 0;

	//DWORD t1 = timeGetTime();
	for(k=0; k<level; k++)
	{
		double scale = pow(zs, k);
		int swd = wd / scale;
		int sht = ht / scale;
		printf("%d %d \n", swd, sht);
		IplImage* pSmall = cvCreateImage( cvSize(swd, sht), pFrame->depth, pFrame->nChannels);
		cvResize(pFrame, pSmall);
		//cvSaveImage("d:\\output\\size.jpg", pSmall);

		//convert from bgr to hsv
		IplImage* hsvFrame = bgr2hsv(pSmall);
		IplImage* h = cvCreateImage( cvGetSize(hsvFrame), IPL_DEPTH_32F, 1 );
		IplImage* s = cvCreateImage( cvGetSize(hsvFrame), IPL_DEPTH_32F, 1 );
		IplImage* v = cvCreateImage( cvGetSize(hsvFrame), IPL_DEPTH_32F, 1 );		
		cvCvtPixToPlane( hsvFrame, h, s, v, NULL );

		//search the resized image
		int nstep = max(2, objWd*0.5/scale);
		//int nstep = 4;
		for(j=0; j<(sht-objHt); j+=nstep)
			for(i=0; i<(swd-objWd); i+=nstep)
			{
				//printf("%d %d . ", j, i);
				int l,r,t,b;
				l = i;
				r = min( l+objWd, swd-1);
				t = j;
				b = min( t+objHt, sht-1);

				//calculate histogram of each window
				histogram hist;
				hist.n = NH*NS + NV;
				memset(hist.histo, 0, sizeof(float)*hist.n);
				for(int kj=t; kj<b; kj++)
					for(int ki=l; ki<r; ki++)
					{
						int bin = histo_bin( 
							((float*)(h->imageData + h->widthStep*kj) )[ki],
							((float*)(s->imageData + s->widthStep*kj) )[ki],
							((float*)(v->imageData + v->widthStep*kj) )[ki] );
						hist.histo[bin] += 1;
					}
					normalize_histogram(&hist);

					//calculate the similarity
					double s = 1- histo_dist_sq(pDstHist, &hist);
					if(s>0.85)
					{
						//printf("%d %d %d %d %lf %lf\n", l, r, t, b, scale, s);
						//cvDrawRect(pDisp, cvPoint(l*scale, (t)*scale), cvPoint(r*scale, (b)*scale), CV_RGB(255,255,0), 2);
						//vecRect.push_back( cvRect(l*scale, t*scale, (r-l)*scale, (b-t)*scale) );
						if(nwindow<1000)
						{
							pRect[nwindow].x = l*scale;
							pRect[nwindow].y = t*scale;
							pRect[nwindow].width = (r-l)*scale;
							pRect[nwindow].height = (b-t)*scale;
						}
						nwindow++;
					}
			}
			cvReleaseImage(&pSmall);
			cvReleaseImage(&hsvFrame);
			cvReleaseImage(&h);
			cvReleaseImage(&s);
			cvReleaseImage(&v);
			//printf("level %d \n", k);
	}
	//DWORD t2 = timeGetTime();

	return nwindow;
}


void VideoSearch(IplImage* pFrame, histogram* pDstHist, vector<CvRect>& vecRect )
{     
	//load video frame image and convert into HSV image
	double ratio = 1; //(double)(pFrame->height) / (double)(pFrame->width);
	int ht = pFrame->height;
	int wd = pFrame->width;
	double zs = 1.3;
	int    objWd = 32;
	int    objHt = objWd*ratio;
	int    level = log(wd/(objWd*4.0)) / log(zs);//log2( (wd/32) );
	int    i,j,k;
	int    nwindow = 0;

	//DWORD t1 = timeGetTime();
	for(k=0; k<level; k++)
	{
		double scale = pow(zs, k);
		int swd = wd / scale;
		int sht = ht / scale;
		printf("%d %d \n", swd, sht);
		IplImage* pSmall = cvCreateImage( cvSize(swd, sht), pFrame->depth, pFrame->nChannels);
		cvResize(pFrame, pSmall);
		//cvSaveImage("d:\\output\\size.jpg", pSmall);

		//convert from bgr to hsv
		IplImage* hsvFrame = bgr2hsv(pSmall);
		IplImage* h = cvCreateImage( cvGetSize(hsvFrame), IPL_DEPTH_32F, 1 );
		IplImage* s = cvCreateImage( cvGetSize(hsvFrame), IPL_DEPTH_32F, 1 );
		IplImage* v = cvCreateImage( cvGetSize(hsvFrame), IPL_DEPTH_32F, 1 );		
		cvCvtPixToPlane( hsvFrame, h, s, v, NULL );

		//search the resized image
		int nstep = max(2, objWd*0.5/scale);
		//int nstep = 4;
		for(j=0; j<(sht-objHt); j+=nstep)
			for(i=0; i<(swd-objWd); i+=nstep)
			{
				//printf("%d %d . ", j, i);
				int l,r,t,b;
				l = i;
				r = min( l+objWd, swd-1);
				t = j;
				b = min( t+objHt, sht-1);

				//calculate histogram of each window
				histogram hist;
				hist.n = NH*NS + NV;
				memset(hist.histo, 0, sizeof(float)*hist.n);
				for(int kj=t; kj<b; kj++)
					for(int ki=l; ki<r; ki++)
					{
						int bin = histo_bin( 
							((float*)(h->imageData + h->widthStep*kj) )[ki],
							((float*)(s->imageData + s->widthStep*kj) )[ki],
							((float*)(v->imageData + v->widthStep*kj) )[ki] );
						hist.histo[bin] += 1;
					}
					normalize_histogram(&hist);

					//calculate the similarity
					double s = 1- histo_dist_sq(pDstHist, &hist);
					if(s>0.85)
					{
						//printf("%d %d %d %d %lf %lf\n", l, r, t, b, scale, s);
						//cvDrawRect(pDisp, cvPoint(l*scale, (t)*scale), cvPoint(r*scale, (b)*scale), CV_RGB(255,255,0), 2);
						vecRect.push_back( cvRect(l*scale, t*scale, (r-l)*scale, (b-t)*scale) );
						nwindow++;
					}
			}
		cvReleaseImage(&pSmall);
		cvReleaseImage(&hsvFrame);
		cvReleaseImage(&h);
		cvReleaseImage(&s);
		cvReleaseImage(&v);
		//printf("level %d \n", k);
	}
	//DWORD t2 = timeGetTime();
}

void VideoSearchUsingIntegralHist(IplImage* pFrame, histogram* pDstHist, vector<CvRect>& vecRect )
{
	int i,j,k;
	double ratio = 1; //(double)(pFrame->height) / (double)(pFrame->width);
	int ht = pFrame->height;
	int wd = pFrame->width;
	double zs = 1.3;
	int    objWd = 32;
	int    objHt = objWd*ratio;
	int    level = log(wd/(objWd*4.0)) / log(zs);//log2( (wd/32) );

	//3. generate integral histogram
	//convert from bgr to hsv
	IplImage* hsvFrame = bgr2hsv(pFrame);
	IplImage* h = cvCreateImage( cvGetSize(hsvFrame), IPL_DEPTH_32F, 1 );
	IplImage* s = cvCreateImage( cvGetSize(hsvFrame), IPL_DEPTH_32F, 1 );
	IplImage* v = cvCreateImage( cvGetSize(hsvFrame), IPL_DEPTH_32F, 1 );		
	cvCvtPixToPlane( hsvFrame, h, s, v, NULL );
	int  scanWdSingle = h->widthStep;

	IntegralHistogram* pIntegralHistImage = (IntegralHistogram*)malloc( (ht+1)*(wd+1)*sizeof(IntegralHistogram) ) ;
	memset(pIntegralHistImage, 0, (ht+1)*(wd+1)*sizeof(IntegralHistogram));

	//clock_t start, finish;
	//double runtime;
	//start  = clock(); 
	GenerateColorIntegralHist( (h->imageData), 
		(s->imageData),
		(v->imageData), 
		ht, wd, scanWdSingle, pIntegralHistImage);
	//finish = clock();		
	//runtime = (double)(finish-start)/CLOCKS_PER_SEC*1000; 
	//printf("\n time for calculate integral histogram: %lf \n", runtime);

	//4. search all scales
	//start  = clock(); 
	for(k=0; k<level; k++)
	{
		//printf("level: %d \n", k);
		double scale = pow(zs, k);
		int swd = wd / scale;
		int sht = ht / scale;
		printf("%d %d \n", swd, sht);

		//search the resized image
		int nstep = max(2, objWd*0.5/scale);
		//int nstep = 4;
		for(j=0; j<(sht-objHt); j+=nstep)
			for(i=0; i<(swd-objWd); i+=nstep)
			{
				//printf("%d %d . ", j, i);
				int l,r,t,b;
				l = i*scale;
				r = min( l+objWd, swd-1)*scale;
				t = j*scale;
				b = min( t+objHt, sht-1)*scale;

				//calculate histogram of each window
				histogram hist;
				hist.n = NH*NS + NV;
				memset(hist.histo, 0, sizeof(float)*hist.n);

				int rbIndex, ltIndex, lbIndex, rtIndex;	
				rbIndex = (b+1)*(wd+1) + r+1;
				ltIndex = t*(wd+1) + l;
				lbIndex = (b+1)*(wd+1) + l;
				rtIndex = t*(wd+1) + r+1;	
				for(int ki=0; ki<hist.n; ki++)
				{		
					hist.histo[ki] = pIntegralHistImage[rbIndex].hist[ki]+pIntegralHistImage[ltIndex].hist[ki]
					-pIntegralHistImage[lbIndex].hist[ki]-pIntegralHistImage[rtIndex].hist[ki];
				}
				normalize_histogram(&hist);

				//calculate the similarity
				double s = 1- histo_dist_sq(pDstHist, &hist);
				if(s>0.85)
				{
					//printf("%d %d %d %d %lf %lf\n", l, r, t, b, scale, s);
					//cvDrawRect(pDisp, cvPoint(l, t), cvPoint(r, b), CV_RGB(255,255,0), 2);
					//nwindow++;
					vecRect.push_back( cvRect(l,t, r-l, b-t) );
				}
			}
	}

}


