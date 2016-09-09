
#include "stdio.h"


#include "CommonFuncs1.h"

//#include "Corelib/commondata.h"
#include "Corelib/CommonFuncs.h"
#include "Corelib/ImageSegment.h"
#include "Corelib/ImageFunc.h"

#include<vector>
using namespace std;


void  VerticalHaarSeg(IplImage* pSrc, IplImage* pDst)
{
	int i,j;

	//segmentation
	int ht = pSrc->height;
	int wd = pSrc->width;

	unsigned char* pSeg = (unsigned char*)malloc(ht*wd);
	//memcpy(pSeg, pbuffer, ht*wd);
	for( j=0; j<ht; j++)
		for( i=0; i<wd; i++)
		{
			pSeg[j*wd+i] = (unsigned char)(pSrc->imageData[j*pSrc->widthStep+i]);
		}
	
	ImageSegOnVerticalHaar(pSeg, ht, wd);
	
	for( j=0; j<ht; j++)
		for(int i=0; i<wd; i++)
		{
			pDst->imageData[j*pDst->widthStep+i] = pSeg[j*wd+i];				
		}

	free(pSeg);
}


double BilinearPixel(double ix, double iy, IplImage* pImage)
{
	int i;
	double   dx,dy;
	int      tx,ty;
	double   w[4];
	double   ng[4];
	double   rg;
	int ht,wd,scanwd;

	ht = pImage->height;
	wd = pImage->width;
	scanwd = pImage->widthStep;

	ix = max(0, min(wd-2, ix));
	iy = max(0, min(ht-2, iy));

	tx = ix;
	ty = iy;
	dx = ix - int(ix);
	dy = iy - int(iy);

	w[0] = (1-dx)*(1-dy);
	w[1] = dx*(1-dy);
	w[2] = (1-dx)*dy;
	w[3] = dx*dy;

	ng[0] = (unsigned char)(pImage->imageData[ty*scanwd+tx]);
	ng[1] = (unsigned char)(pImage->imageData[ty*scanwd+tx+1]);
	ng[2] = (unsigned char)(pImage->imageData[(ty+1)*scanwd+tx]); 
	ng[3] = (unsigned char)(pImage->imageData[(ty+1)*scanwd+tx+1]);

	rg = 0;
	for(i=0; i<4; i++)
		rg += w[i]*ng[i];

	return rg;
}


/*  convert IplImage to float image
*/
void   IplImageToFloatImage(IplImage* pImage, float** pBuffer, int* ht, int* wd)
{

	if(pImage==NULL)
		return;

	*ht = pImage->height;
	*wd = pImage->width;
	*pBuffer = (float*)malloc( (*ht)*(*wd)*sizeof(float) );

	//8 bits
	if(pImage->nChannels == 1)
	{
		int scanwd = pImage->widthStep;
		for(int j=0; j<(*ht); j++)
			for(int i=0; i<(*wd); i++)
			{
				(*pBuffer)[j*(*wd)+i] = (unsigned char)( pImage->imageData[j*scanwd+i] );
			}
	}

	//24 bits
	if(pImage->nChannels == 3)
	{
		IplImage* pGray = cvCreateImage( cvGetSize(pImage), 8, 1);
		int scanwd = pGray->widthStep;
		cvCvtColor(pImage, pGray, CV_BGR2GRAY);		

		for(int j=0; j<(*ht); j++)
			for(int i=0; i<(*wd); i++)
			{
				(*pBuffer)[j*(*wd)+i] = (unsigned char)(pGray->imageData[j*scanwd+i]);
			}
			cvReleaseImage(&pGray);
	}
}



/*  convert IplImage to gray image
*/
void   IplImageToGrayImage(IplImage* pImage, unsigned char** pBuffer, int* ht, int* wd)
{
	if(pImage==NULL)
		return;

	*ht = pImage->height;
	*wd = pImage->width;
	*pBuffer = (unsigned char*)malloc( (*ht)*(*wd) );

	//8 bits
	if(pImage->nChannels == 1)
	{
		int scanwd = pImage->widthStep;
		for(int j=0; j<(*ht); j++)
			for(int i=0; i<(*wd); i++)
			{
				(*pBuffer)[j*(*wd)+i] = (unsigned char)( pImage->imageData[j*scanwd+i] );
			}
	}

	//24 bits
	if(pImage->nChannels == 3)
	{
		IplImage* pGray = cvCreateImage( cvGetSize(pImage), 8, 1);
		int scanwd = pGray->widthStep;
		cvCvtColor(pImage, pGray, CV_BGR2GRAY);		
		
		for(int j=0; j<(*ht); j++)
			for(int i=0; i<(*wd); i++)
			{
				(*pBuffer)[j*(*wd)+i] = (unsigned char)(pGray->imageData[j*scanwd+i]);
			}
		cvReleaseImage(&pGray);
	}
}


/*  convert IplImage to gray image
*/
void   IplImageToColorImage(IplImage* pImage, unsigned char** r, unsigned char** g, unsigned char** b, int* ht, int* wd)
{
	if(pImage==NULL)
		return;

	*ht = pImage->height;
	*wd = pImage->width;
	*r = (unsigned char*)malloc( (*ht)*(*wd) );
	*g = (unsigned char*)malloc( (*ht)*(*wd) );
	*b = (unsigned char*)malloc( (*ht)*(*wd) );

	//8 bits
	if(pImage->nChannels == 1)
	{
		int scanwd = pImage->widthStep;
		for(int j=0; j<(*ht); j++)
			for(int i=0; i<(*wd); i++)
			{
				(*r)[j*(*wd)+i] = (unsigned char)( pImage->imageData[j*scanwd+i] );
				(*g)[j*(*wd)+i] = (unsigned char)( pImage->imageData[j*scanwd+i] );
				(*b)[j*(*wd)+i] = (unsigned char)( pImage->imageData[j*scanwd+i] );
			}
	}

	//24 bits
	if(pImage->nChannels == 3)
	{
		int scanwd = pImage->widthStep;
		for(int j=0; j<(*ht); j++)
			for(int i=0; i<(*wd); i++)
			{
				(*r)[j*(*wd)+i] = (unsigned char)( pImage->imageData[j*scanwd+i*3] );
				(*g)[j*(*wd)+i] = (unsigned char)( pImage->imageData[j*scanwd+i*3+1] );
				(*b)[j*(*wd)+i] = (unsigned char)( pImage->imageData[j*scanwd+i*3+2] );
			}		
	}
}

void GenerateOpticalFlow(IplImage* pFirst, IplImage* pSecond, char* filepath)
{
	IplImage *curPYR = 0;
	IplImage *nextPYR = 0;
	CvPoint tp1;
	MyPointF pt;
	//int len=3;
	CvScalar color;
	int ht,wd,scanwd;
	int Level = 2;
	int k = 0;
	int i,j;
	double* px;
	double* py;		
	int     npt;	
	int     index[3];
	double  tx,ty;
	MyPointF tpt[3];
	double len[3];
	double weight[3];
	MyPointF tp;
	double alllen;			

	ht = pFirst->height;
	wd = pFirst->width;
	scanwd = pFirst->widthStep;

	curPYR  = cvCreateImage( cvSize(wd, ht), 8, 1 );
	nextPYR = cvCreateImage( cvSize(wd, ht), 8, 1 );

	//1. find good features (harris feature points)
	IplImage* eig = cvCreateImage( cvGetSize(pFirst), 32, 1 );
	IplImage* temp = cvCreateImage( cvGetSize(pFirst), 32, 1 );
	double quality = 0.01;
	double min_distance = 10;
	int win_size = 10;
	CvPoint2D32f* points = (CvPoint2D32f*)malloc(MAX_COUNT*sizeof(CvPoint2D32f));
	CvPoint2D32f* npoints = (CvPoint2D32f*)malloc(MAX_COUNT*sizeof(CvPoint2D32f));        
	int count = MAX_COUNT;
	cvGoodFeaturesToTrack( pFirst, eig, temp, points, &count,
		quality, min_distance, 0, 3, 0, 0.04 );
	//cvFindCornerSubPix( pImage, points, count,
	//	cvSize(win_size,win_size), cvSize(-1,-1),
	//	cvTermCriteria(CV_TERMCRIT_ITER|CV_TERMCRIT_EPS,20,0.03));
	//save the feature points from bench image 	
	cvReleaseImage( &eig );
	cvReleaseImage( &temp );

	//2. calculate optical flow
	char* status = (char*)malloc(count);
	cvCalcOpticalFlowPyrLK( pFirst, pSecond, curPYR, nextPYR, 
		points, npoints, 
		count, 
		cvSize(10,10), 
		Level, 
		status, 0, 
		cvTermCriteria(CV_TERMCRIT_ITER|CV_TERMCRIT_EPS,20,0.03),
		0);
	free(status);

	//calculate the mean translation		
	tx = 0;
	ty = 0;
	for(k=0; k<count; k++)
	{
		tx += (points[k].x - npoints[k].x);
		ty += (points[k].y - npoints[k].y);
	}
	tx /= (double)(count);
	ty /= (double)(count);

	//save the sparse optical flow 
	//char sopticalFlow[256];
	//strcpy(sopticalFlow, m_binPath);
	//strcat(sopticalFlow, "opticalFlow.txt");

	FILE* fp = fopen(filepath, "w");
	fprintf(fp, "%d \n", count);
	for(k=0; k<count; k++)
	{
		if( npoints[k].x<-3 || npoints[k].x>wd || npoints[k].y<-3 || npoints[k].y>ht  )
		{
			continue;
		}

		fprintf(fp, "%lf %lf %lf %lf \n", points[k].x, points[k].y, npoints[k].x,npoints[k].y); 
	}
	fclose(fp);

	cvReleaseImage(&curPYR);
	cvReleaseImage(&nextPYR);
	free(points);
	free(npoints);
}


void ConnectedComonent(IplImage* pImage, MyRect* pRect, int* nRect)
{
	int ht = pImage->height;
	int wd = pImage->width;
	int scanwd = pImage->widthStep;
	unsigned char* mask = (unsigned char*)malloc(ht*wd);
	for(int j=0; j<ht; j++)
		for(int i=0; i<wd; i++)
		{
			mask[j*wd+i] = pImage->imageData[j*scanwd+i];
		}
	*nRect = 0;
	SegLabel(mask, ht, wd, pRect, nRect);
	free(mask);
}


void PlateImageSeg(IplImage* pImage, unsigned char* mask)
{
	IplImage* pInput = cvCloneImage(pImage);

	cvSobel(pInput, pInput,1,0);
	cvThreshold(pInput, pInput, 30, 255, CV_THRESH_BINARY);
	cvSmooth(pInput, pInput, CV_MEDIAN);	

	//2. Morphological filtering according to the characteristics of plate
	IplConvKernel* element = NULL;
	element = cvCreateStructuringElementEx( 15, 1, 7, 0, CV_SHAPE_RECT, 0 );		
	cvDilate(pInput, pInput, element);
	cvReleaseStructuringElement(&element);   
	//belement = cvCreateStructuringElementEx( 7, 1, 1, 1, CV_SHAPE_RECT, 0 );
	cvErode(pInput, pInput, element);
	cvReleaseStructuringElement(&element);					

	//3. connected components
	int ht = pInput->height;
	int wd = pInput->width;
	int scanwd = pInput->widthStep;
	for(int j=0; j<ht; j++)
		for(int i=0; i<wd; i++)
		{
			mask[j*wd+i] = pInput->imageData[j*scanwd+i];
		}

	cvSaveImage("d:\\mask.jpg", pInput);

	cvReleaseImage(&pInput);
}

void VerticalEdgeDetect(IplImage* pSrc, IplImage* pDst)
{
	int ht = pSrc->height;
	int wd = pSrc->width;
	int scanWd = pSrc->widthStep;

	//1. sobel gradient
	cvSmooth(pSrc, pSrc, CV_GAUSSIAN);
	IplImage* pSobel = cvCloneImage(pSrc);
	cvSobel(pSrc, pSobel, 1, 0);
	//cvSmooth(pSobel, pSobel, CV_MEDIAN);

	//2. vertical edge detection		
	IplImage* pEdge = cvCloneImage(pSobel);
	memset(pEdge->imageData, 0, ht*scanWd);
	for(int j=0; j<ht; j++)
		for(int i=2; i<wd-2; i++)
		{
			int c = (unsigned char)( pSobel->imageData[j*scanWd+i] );
			int l = (unsigned char)( pSobel->imageData[j*scanWd+i-1] );
			int r = (unsigned char)( pSobel->imageData[j*scanWd+i+1] );

			if( c>30 )
			{					
				if( c>l && c>=r)
				{
					pDst->imageData[j*scanWd+i] = 255;
				}
				else
				{
					pDst->imageData[j*scanWd+i] = 0;
				}					
			}
			else
			{
				pDst->imageData[j*scanWd+i] = 0;
			}
		}
}

void IplImageSplit(IplImage* pSrc, IplImage* pR, IplImage* pG, IplImage* pB)
{
	int ht = pSrc->height;
	int wd = pSrc->width;
	int scanWd = pSrc->widthStep;

	if(pSrc->nChannels!=3)
		return;

	for(int j=0; j<ht; j++)
		for(int i=0; i<wd; i++)
		{
			pR->imageData[j*wd+i] = pSrc->imageData[j*scanWd+i*3];
			pG->imageData[j*wd+i] = pSrc->imageData[j*scanWd+i*3+1];
			pB->imageData[j*wd+i] = pSrc->imageData[j*scanWd+i*3+2];
		}
}



double SegDistance(stLINE l1, stLINE l2)
{
	double   len = 0;
	double   len1,len2;
	MyPointF p1,p2;
	p1.x = (l1.v1.x + l1.v2.x)*0.5;
	p1.y = (l1.v1.y + l1.v2.y)*0.5;
	p2.x = (l2.v1.x + l2.v2.x)*0.5;
	p2.y = (l2.v1.y + l2.v2.y)*0.5;
	len = P2PDistance(p1, p2);

	len1 = P2PDistance( l1.v1, l1.v2 );
	len2 = P2PDistance( l2.v1, l2.v2 );

	double d = len - (len1+len2)*0.5;

	return d;
}

void FitLine(MyPointF* pPts, int count, double* a, double* b, double* c)
{

	CvPoint* points = (CvPoint*)malloc( count * sizeof(points[0]));
	CvMat pointMat = cvMat( 1, count, CV_32SC2, points );
	float line[4];

	for(int i=0; i<count; i++)
	{
		points[i].x = pPts[i].x;
		points[i].y = pPts[i].y;
	}

	//
	cvFitLine( &pointMat, CV_DIST_L1, 1, 0.001, 0.001, line );

	double vx = line[0];
	double vy = line[1];
	double x0 = line[2];
	double y0 = line[3];
	*a = vy;
	*b = -vx;
	*c = -vy*x0 + vx*y0;

	free(points);
}

stLINE MergeTwoParallelLines(stLINE l1, stLINE l2)
{
	stLINE   m;
	
	//
	//m.rou = (l1.rou + l2.rou)*0.5;
	//m.sita = (l1.sita + l2.sita)*0.5;
    MyPointF pts[4];
    pts[0] = l1.v1;
	pts[1] = l1.v2;
	pts[2] = l2.v1;
	pts[3] = l2.v2;
	double a,b,c;
	FitLine(pts, 4, &a, &b, &c);
	m.a = a; m.b = b; m.c = c;
		
	//calculate the end vertex
	double minx,maxx;
	double miny,maxy;
	
	minx = min( min(l1.v1.x, l1.v2.x),  min(l2.v1.x, l2.v2.x) );
	maxx = max( max(l1.v1.x, l1.v2.x),  max(l2.v1.x, l2.v2.x) );	
	miny = -(a*minx+c)/b;
	maxy = -(a*maxx+c)/b;
	/*
	double angle = m.sita/180.0*PI;
	miny = (m.rou - minx*cos(angle)) / sin(angle);
	maxy = (m.rou - maxx*cos(angle)) / sin(angle);
	*/
	m.v1.x = minx;
	m.v1.y = miny;
	m.v2.x = maxx;
	m.v2.y = maxy;
    CalculateLinePolarParams(&m);


	return m;
}

/*
void MergeParallelLines(vector<stLINE>& lines)
{
	int i,j;
	int num;

	vector<stLINE> mergeLines;
	int* pmask = (int*)malloc( lines.size()*sizeof(int) );
	memset(pmask, 0, sizeof(int)*lines.size());

	num = lines.size();
	for(j=0; j<num; j++)
	{
		if(pmask[j]>0)
			continue;

		stLINE mline = lines[j];
		for(i=0; i<num; i++)
		{
			if(pmask[i]>0)
				continue;

			//if the two lines are parallel and close enough
			if(  fabs(mline.rou-lines[i].rou)<5 && fabs(mline.sita-lines[i].sita)<5 ) 
			{
				int nIsMerge =  MergeTwoParallelLines( mline,  lines[i], mline);
				if(nIsMerge)
				{  
					pmask[i]=1;
				}
			}			
		}
		mergeLines.push_back(mline);
		pmask[j] = 1;
	}

	//overwrite
	lines = mergeLines;
}
*/

void MergeParallelLines(vector<stLINE>& lines)
{
	int i,j;
	int num;

	vector<stLINE> mergeLines;
	int* pmask = (int*)malloc( lines.size()*sizeof(int) );
	memset(pmask, 0, sizeof(int)*lines.size());

	num = lines.size();
	for(j=0; j<num; j++)
	{
		if(pmask[j]>0)
			continue;

		stLINE mline = lines[j];		
		//find the closet line
		while(1) 
		{
			double minD = 10000000;
			int    index = 0;
			for(i=0; i<num; i++)
			{
				if(i==j)
					continue;
				if(pmask[i]>0)
					continue;

				if( fabs(mline.rou-lines[i].rou)<5 && fabs(mline.sita-lines[i].sita)<5 ) 
				{
					double d = SegDistance(mline, lines[i]);
					if(d<minD)
					{
						minD = d;
						index = i;
					}
				}			
			}
			//if the two lines are parallel and close enough
			if(minD<20)
			{
				mline = MergeTwoParallelLines(mline, lines[index]);
				pmask[index] = 1;				
			}
			else
			{
				break;
			}
		} 	

		mergeLines.push_back(mline);
		pmask[j] = 1;
	}

	//overwrite
	lines = mergeLines;
}