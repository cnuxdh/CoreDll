
#include "stdio.h"
#include "string.h"
#include "MosaicFuncsV1.h"
//#include "image.h"
#include "CommonFuncs.h"


#include "../../Corelib/commonfile.h"
#include "../../Corelib/matching.h"
#include "../../Corelib/ImageFunc.h"
#include "../../Corelib/commondata.h"
#include "../../Corelib/Triangles.h"
#include "../../Corelib/CommonFuncs.h"
#include "../../Corelib/image.h"
#include "../../Corelib/Matrix.h"
#include "../../Corelib/FitObject.h"


#include "defs.h"
#include "vim_imgs.h"
#include "vim_sift.h"
#include "ransac.h"

#include <vector>
using namespace std;


#define NN_SQ_DIST_RATIO   0.8
#define MAX_FEAT_NUM 10000


/*
int FundamentalRefine(MyPointF* pLeftPt, MyPointF* pRightPt, int npt, CvMat** fundamental_matrix)
{
	CvMat* points1;
	CvMat* points2;
	CvMat* status;
	
	int i;

	points1 = cvCreateMat(1,npt,CV_32FC2);
	points2 = cvCreateMat(1,npt,CV_32FC2);
	status  = cvCreateMat(1,npt,CV_8UC1);

	//Fill the points here
	for( i = 0; i < npt; i++ )
	{
		points1->data.fl[i*2]   = pLeftPt[i].x;  //These are points such as found
		points1->data.fl[i*2+1] = pLeftPt[i].y;  // on the chessboard calibration
		points2->data.fl[i*2]   = pRightPt[i].x;  // pattern.
		points2->data.fl[i*2+1] = pRightPt[i].y;
	}
	
	*fundamental_matrix = cvCreateMat(3,3,CV_32FC1);

	int fm_count = cvFindFundamentalMat( points1, points2,
					*fundamental_matrix,
					CV_FM_RANSAC, 1.0, 0.99, status );

	int nInlier=0;

	for(i=0; i<npt; i++)
	{
		if(status->data.ptr[i]>0)
		{
			pLeftPt[nInlier].x = points1->data.fl[i*2];
			pLeftPt[nInlier].y = points1->data.fl[i*2+1];

			pRightPt[nInlier].x = points2->data.fl[i*2];
			pRightPt[nInlier].y = points2->data.fl[i*2+1];

			nInlier++;
		}
	}

	cvReleaseMat(&points1);
	cvReleaseMat(&points2);
	//cvReleaseMat(&fundamental_matrix);
	cvReleaseMat(&status);
	return nInlier;
}



//去掉冗余的匹配点对，即一对多，或者多对一
int DeleteRedundantPairs(MyPointF* pt1, MyPointF* pt2, int np)
{		
	int i,j,k;
	MyPointF* tp1 = NULL;
	MyPointF* tp2 = NULL;
	int tnp;
	int* mask = NULL;

	mask = (int*)malloc(sizeof(int)*np);
	memset(mask, 0, sizeof(int)*np);

	tp1 = (MyPointF*)malloc( np*sizeof(MyPointF) );
	tp2 = (MyPointF*)malloc( np*sizeof(MyPointF) );

	for( i=0; i<np; i++ )
	{
		for(j=i+1; j<np; j++)
		{
			if( pt1[i].x==pt1[j].x 
				&& pt1[i].y==pt1[j].y )
			{
				mask[i] = 1;
				mask[j] = 1;
			}
		}
	}
	for( i=0; i<np; i++ )
	{
		for(j=i+1; j<np; j++)
		{
			if( pt2[i].x==pt2[j].x 
				&& pt2[i].y==pt2[j].y )
			{
				mask[i] = 1;
				mask[j] = 1;
			}
		}
	}

	tnp = 0;
	for(i=0; i<np; i++)
	{
		if(mask[i]==0)
		{
			tp1[tnp] = pt1[i];
			tp2[tnp] = pt2[i];
			tnp ++;
		}
	}
	for(i=0; i<tnp; i++)
	{
		pt1[i] = tp1[i];
		pt2[i] = tp2[i];		
	}

	free(mask);
	free(tp1);
	free(tp2);
	return tnp;
}


//删除特征点数组中的位置相同的点
int DeleteSamePoints(MyPointF* feat, int n)
{
	MyPointF* temp;
	int i,j;
    int *mask;
	int np;
	
	temp = new MyPointF[n];
	
	mask = (int*)malloc(n*sizeof(int));
	memset(mask, 0, sizeof(int)*n);
	
	for(i=0; i<n; i++)
	{
		for(j=i+1; j<n; j++)
		{
			if( mask[j]==1 )
				continue;
			
			if(   ( (feat+i)->x==(feat+j)->x) 
				&& ( (feat+i)->y==(feat+j)->y) )
			{
				mask[j] = 1;
			}
		}
	}
	
	np = 0;
	for(i=0; i<n; i++)
	{
		if(mask[i]==0)
		{
			temp[np] = *(feat+i);			
			np ++;
		}
	}
	for(i=0; i<np; i++)
	{
		*(feat+i) = temp[i];		
	}
	
	free(mask);
	delete[] temp;
	return np;
}


//删除特征点数组中的位置相同的点
int DeleteSamePoints(struct feature* feat, int n)
{
	struct feature* temp;
	int i,j;
    int *mask;
	int np;

	temp = new struct feature[n];

	mask = (int*)malloc(n*sizeof(int));
	memset(mask, 0, sizeof(int)*n);

	for(i=0; i<n; i++)
	{
		for(j=i+1; j<n; j++)
		{
			if( mask[j]==1 )
				continue;

			if(   ( (feat+i)->x==(feat+j)->x) 
				&& ( (feat+i)->y==(feat+j)->y) )
			{
				mask[j] = 1;
			}
		}
	}
		
	np = 0;
	for(i=0; i<n; i++)
	{
		if(mask[i]==0)
		{
			temp[np] = *(feat+i);			
			np ++;
		}
	}
	for(i=0; i<np; i++)
	{
		*(feat+i) = temp[i];		
	}
	
	free(mask);
	delete[] temp;
	return np;
}


//加载sift特征数据文件
void LoadSiftFile(char* filename, struct feature** feat, int* n)
{
	FILE* fp = NULL;
	*n = 0;
	
	fp = fopen(filename, "rb");
	fread(n, 1, sizeof(int), fp);
	*feat = (struct feature*)malloc(*n*sizeof(struct feature));
	fread(*feat, *n, sizeof(struct feature), fp);
	fclose(fp);	
}
*/

//////////////////////////////////////////////////////////////////////////

/*
Combines two images by scacking one on top of the other
img1 - top image
img2 - bottom image
return : 
Returns the image resulting from stacking \a img1 on top if \a img2
*/
IplImage* stack_imgs( IplImage* img1, IplImage* img2 )
{
	IplImage* stacked = cvCreateImage( cvSize( MAX(img1->width, img2->width), \
		img1->height + img2->height ), \
		IPL_DEPTH_8U, 3 );
	cvZero( stacked );
	cvSetImageROI( stacked, cvRect( 0, 0, img1->width, img1->height ) );
	cvAdd( img1, stacked, stacked, NULL );
	cvSetImageROI( stacked, cvRect(0, img1->height, img2->width, img2->height) );
	cvAdd( img2, stacked, stacked, NULL );
	cvResetImageROI( stacked );	
	return stacked;
}

IplImage* MosaicTwoImagesPolynomial(IplImage* srcImage, IplImage* dstImage, float* p)
{

	int i,j,k;
	IplImage* mosaicImage = NULL;
	CvScalar s;
	int tx,ty;
    int cx[4];
	int cy[4];
	int ht,wd;

	int minX, minY, maxX, maxY;

	ht = srcImage->height;
	wd = srcImage->width;

	minX = 0;
	minY = 0;
	maxX = wd;
	maxY = ht;

	cx[0] = 0;  cy[0] = 0;
	cx[1] = wd; cy[1] = 0;
	cx[2] = 0;  cy[2] = ht;
	cx[3] = wd; cy[3] = ht;

	//calculate the polygon of src image in the dst image
	for(i=0; i<4; i++)
	{
		tx = POLYNOMIAL_X(cx[i], cy[i], p);
		ty = POLYNOMIAL_Y(cx[i], cy[i], p);		
		minX = min(minX, tx);
		minY = min(minY, ty);
		maxX = max(maxX, tx);
		maxY = max(maxY, ty);
	}

	mosaicImage = cvCreateImage( cvSize(maxX-minX, maxY-minY), IPL_DEPTH_8U, 3);

	for(j=0; j<dstImage->height; j++)
		for(i=0; i<dstImage->width; i++)
		{
			s=cvGet2D(dstImage, j, i);
			cvSet2D(mosaicImage, j-minY, i-minX, s);
		}    

	for(j=0; j<srcImage->height; j++)
		for(i=0; i<srcImage->width; i++)
		{
			s=cvGet2D(srcImage, j, i); // get the (i,j) pixel value
			tx = POLYNOMIAL_X(i,j,p);	
			ty = POLYNOMIAL_Y(i,j,p);	
           
			if( tx>=(maxX-minX) || tx<0 
				|| ty>=(maxY-minY) || ty<0 )
				continue;

			cvSet2D(mosaicImage, ty-minY, tx-minX, s);	
			
			/*
			printf("B=%f, G=%f, R=%f ",s.val[0],s.val[1],s.val[2]);
			s.val[0]=111;
			s.val[1]=111;
			s.val[2]=111;
			cvSet2D(srcImage, i, j, s);//set the (i,j) pixel value
			*/
		}	
	return mosaicImage;

}


//#define AFFINE_X(x,y,p) p[0]*x+p[1]*y+p[2]
//#define AFFINE_Y(x,y,p) p[3]*x+p[4]*y+p[5]

/* function: mosaic srcImage into dstImage based on affine transform
*/
IplImage* MosaicTwoImages(IplImage* srcImage, IplImage* dstImage, float* p, int* nminx, int* nminy)
{
	int i,j,k;
	IplImage* mosaicImage = NULL;
	CvScalar s;
	int tx,ty;
    int cx[4];
	int cy[4];
	int ht,wd;
	float ip[6];

	int minX, minY, maxX, maxY;
	int cminx,cminy,cmaxx,cmaxy;

	ht = srcImage->height;
	wd = srcImage->width;

	minX = 0;
	minY = 0;
	maxX = dstImage->width;
	maxY = dstImage->height;

	cx[0] = 0;  cy[0] = 0;
	cx[1] = wd; cy[1] = 0;
	cx[2] = 0;  cy[2] = ht;
	cx[3] = wd; cy[3] = ht;

	//calculate the polygon of src image in the dst image
	for(i=0; i<4; i++)
	{
		tx = AFFINE_X(cx[i], cy[i], p);
		ty = AFFINE_Y(cx[i], cy[i], p);		
		minX = min(minX, tx);
		minY = min(minY, ty);
		maxX = max(maxX, tx);
		maxY = max(maxY, ty);

		cx[i] = tx;
		cy[i] = ty;
	}

	mosaicImage = cvCreateImage( cvSize(maxX-minX, maxY-minY), IPL_DEPTH_8U, 3);

	for(j=0; j<dstImage->height; j++)
		for(i=0; i<dstImage->width; i++)
		{
			s=cvGet2D(dstImage, j, i);
			cvSet2D(mosaicImage, j-minY, i-minX, s);
		}    
		
	/*
	for(j=0; j<srcImage->height; j++)
		for(i=0; i<srcImage->width; i++)
		{
			s=cvGet2D(srcImage, j, i); // get the (i,j) pixel value
			tx = AFFINE_X(i,j,p);	
			ty = AFFINE_Y(i,j,p);		
			
			if( (ty-minY)>=0 && (ty-minY)<mosaicImage->height
				&& (tx-minX)>=0 && (tx-minX)<mosaicImage->width )
				cvSet2D(mosaicImage, ty-minY, tx-minX, s);	
		}
	*/

	cminx = 10000;
	cmaxx = 0;
	cminy = 10000;
	cmaxy = 0;
	for(i=0; i<4; i++)
	{
		if(cminx>cx[i]) cminx=cx[i];
		if(cmaxx<cx[i]) cmaxx=cx[i];
		if(cminy>cy[i]) cminy=cy[i];
		if(cmaxy<cy[i]) cmaxy=cy[i];
	}

	for(i=0; i<6; i++)
		ip[i] = p[i];

	InvertAffine(ip);    

	for(j=cminy; j<cmaxy; j++)
		for(i=cminx; i<cmaxx; i++)
		{
			int ix,iy;
			ix = AFFINE_X(i,j,ip);
			iy = AFFINE_Y(i,j,ip);
			if(ix>0 && ix<wd && iy>0 && iy<ht)
			{
				s=cvGet2D(srcImage, iy, ix); // get the (i,j) pixel value				
				
				if( s.val[0]<10 && s.val[1]<10 && s.val[2]<10 )
					continue;

				cvSet2D(mosaicImage, j-minY, i-minX, s);
			}
		}

	*nminx = minX;
	*nminy = minY;
		
	return mosaicImage;
}

/*
  function: mosaic multiple images based on affine transform
*/
void MosaicMultiImages(char* imagePath, double** pTransPara, int nBegin, int nEnd)
{

}

/*  功能：基于sift特征进行配准
    返回：
	    配准后的点的索引值
*/
void FindCorrespondenceIndex(struct feature* feat1, int n1,
							 struct feature* feat2, int n2,
							 int *pi1, int* pi2)
{

	


}



/* 在两张图片的sift特征之间建立对应关系( 并返回配准点对应的索引值 )
   input: 
     feat1,feat2,n1,n2: 左右两张相片的特征点及其数量
   output:
	 p1,p2,np: 匹配点及其数量
	 ip1,ip2:  匹配点对应的索引值
 */
/*
void FindCorrespondenceNew(struct feature* feat1, int n1,
						struct feature* feat2, int n2,
						MyPointF* p1, MyPointF* p2, 
						int *ip1, int* ip2,
						int* np)
{
	struct feature* feat;
	struct feature** nbrs;
	struct kd_node* kd_root;
	double d0, d1;
	int k, i, m = 0;
	CvPoint2D32f pt1, pt2;

	*np = 0;
	
	m = 0;
	//----------Find the correspondence of two type keypoints ---------
	kd_root = kdtree_build( feat2, n2 );

	for( i = 0; i < n1; i++ )
	{
		feat = feat1 + i;
		k = kdtree_bbf_knn( kd_root, feat, 2, &nbrs, 200 );
		if( k == 2 )
		{
			d0 = descr_dist_sq( feat, nbrs[0] );
			d1 = descr_dist_sq( feat, nbrs[1] );
			if( d0 < d1 * NN_SQ_DIST_RATIO )
			{
				ip1[m] = i;
				ip2[m] = nbrs[0]->index;

				pt1.x = feat->x;
				pt1.y = feat->y;
				pt2.x = nbrs[0]->x;
				pt2.y = nbrs[0]->y;

				//pt1 = CvPoint2D32f( cvRound( feat->x ), cvRound( feat->y ) );
				//pt2 = CvPoint2D32f( cvRound( nbrs[0]->x ), cvRound( nbrs[0]->y ) );
				
				//cvSeqPush(img1_seq,&pt1); 
				//cvSeqPush(img2_seq,&pt2);				
				//pt2.y += img1->height;				
				//cvLine( stacked, pt1, pt2, CV_RGB(255,0,255), 1, 8, 0 );								
				//feat1[i].fwd_match = nbrs[0];

				p1[m].x = pt1.x;
				p1[m].y = pt1.y;
				p2[m].x = pt2.x;
				p2[m].y = pt2.y;
				
				m++;
			}
		}

		free( nbrs );
	}

	*np = m;

	kdtree_release( kd_root );
}
*/


/* 在两张图片的sift特征之间建立对应关系
input: 
feat1,feat2,n1,n2: 左右两张相片的特征点及其数量
output:
p1,p2,np: 匹配点及其数量
*/
/*
void FindCorrespondence(struct feature* feat1, int n1,
						struct feature* feat2, int n2,
						struct feature* p1, struct feature* p2, int* np)
{
	struct feature* feat;
	struct feature** nbrs;
	struct kd_node* kd_root;
	double d0, d1;
	int k, i, m = 0;
	//CvPoint2D32f pt1, pt2;

	*np = 0;

	m = 0;
	//----------Find the correspondence of two type keypoints ---------
	kd_root = kdtree_build( feat2, n2 );
	for( i = 0; i < n1; i++ )
	{
		feat = feat1 + i;
		k = kdtree_bbf_knn( kd_root, feat, 2, &nbrs, 200 );
		if( k == 2 )
		{ 
			d0 = descr_dist_sq( feat, nbrs[0] );
			d1 = descr_dist_sq( feat, nbrs[1] );
			if( d0 < d1 * NN_SQ_DIST_RATIO )
			{
				memcpy( p1+m, feat, sizeof(feature));
				memcpy( p2+m, nbrs[0], sizeof(feature));
				m++;
			}
		}
		free( nbrs );
	}
	*np = m;

	kdtree_release( kd_root );
}
*/

/* 在两张图片的sift特征之间建立对应关系
   input: 
     feat1,feat2,n1,n2: 左右两张相片的特征点及其数量
   output:
	 p1,p2,np: 匹配点及其数量
 */
/*
void FindCorrespondence(struct feature* feat1, int n1,
						struct feature* feat2, int n2,
						MyPointF* p1, MyPointF* p2, int* np)
{
	struct feature* feat;
	struct feature** nbrs;
	struct kd_node* kd_root;
	double d0, d1;
	int k, i, m = 0;
	CvPoint2D32f pt1, pt2;

	*np = 0;
	
	m = 0;
	//----------Find the correspondence of two type keypoints ---------
	kd_root = kdtree_build( feat2, n2 );
	for( i = 0; i < n1; i++ )
	{
		feat = feat1 + i;
		k = kdtree_bbf_knn( kd_root, feat, 2, &nbrs, 200 );
		if( k == 2 )
		{ 
			d0 = descr_dist_sq( feat, nbrs[0] );
			d1 = descr_dist_sq( feat, nbrs[1] );
			if( d0 < d1 * NN_SQ_DIST_RATIO )
			{
				pt1.x = feat->x;
				pt1.y = feat->y;
				pt2.x = nbrs[0]->x;
				pt2.y = nbrs[0]->y;	

				p1[m].x = pt1.x;
				p1[m].y = pt1.y;
				p2[m].x = pt2.x;
				p2[m].y = pt2.y;				
				m++;
			}
		}
		free( nbrs );
	}
	*np = m;
	
	kdtree_release( kd_root );
}
*/

/*
void GenerateTriangle(IplImage* img, MyPointF* pt, int np)
{
	int i,j;
	CTINClass* pTin = NULL;		
	
	pTin = new CTINClass("aaa");			
	pTin->BeginAddPoints();
	for(i=0; i<np; i++)
	{
		pTin->AddPoint(pt[i].x, pt[i].y, i);
	}
	pTin->EndAddPoints();
	pTin->FastConstruct();
	pTin->EnableIntersection( FALSE );
	
	//draw triangle
	LONG nTriangleNum = 0 ;		
	TRIANGLE **tris = pTin->SelectTriangles(&nTriangleNum, 0, 0, 0, 0 );	
	CvPoint pts[3];
	for(i=0; i<nTriangleNum; i++)
	{		
		for(int k=0; k<3; k++)
		{		
			pts[k].x = (*tris)->vertex[k]->x;
			pts[k].y = (*tris)->vertex[k]->y;
		}
		cvLine(img, pts[0], pts[1], CV_RGB(255,0,0), 1, 8, 0);
		cvLine(img, pts[1], pts[2], CV_RGB(255,0,0), 1, 8, 0);
		cvLine(img, pts[0], pts[2], CV_RGB(255,0,0), 1, 8, 0);
		*tris++;	
	}	
	delete pTin;
}


IplImage* MosaicUsingTriangle(IplImage* srcImage, IplImage* dstImage,
						 MyPointF* srtPt, MyPointF* dstPt, int np)
{
	int i,j,k;
	int m,n;
	//generate triangle
	CTINClass* pTin = NULL;		
	CvPoint pts[3];
	CvPoint pts2[3];
	float px[3];
	float py[3];
	int index;
	int minX, minY, maxX, maxY;
	int ht,wd;
	int cx[4],cy[4];
	IplImage* mosaicImage;
	CvScalar s;
	int tx,ty;
	float ap[6];
	LONG nTriangleNum = 0 ;	
	TRIANGLE **tris;
	
	pTin = new CTINClass("aaa");			
	pTin->BeginAddPoints();
	for(i=0; i<np; i++)
	{
		pTin->AddPoint(srtPt[i].x, srtPt[i].y, i);
	}
	pTin->EndAddPoints();
	pTin->FastConstruct();
	pTin->EnableIntersection( FALSE );
	
	//get each triangle		
	tris = pTin->SelectTriangles(&nTriangleNum, 0, 0, 0, 0 );	
	printf("\n triangle number: %d \n", nTriangleNum);	
	for(i=0; i<nTriangleNum; i++)
	{		
		for(k=0; k<3; k++)
		{		
			pts[k].x = (*tris)->vertex[k]->x;
			pts[k].y = (*tris)->vertex[k]->y;
			index = (*tris)->vertex[k]->attr;
			printf(" %d ", index);

			pts2[k].x = dstPt[index].x;
			pts2[k].y = dstPt[index].y;
		}
		printf("\n");

		cvLine(srcImage, pts[0], pts[1], CV_RGB(255,0,0), 1, 8, 0);
		cvLine(srcImage, pts[1], pts[2], CV_RGB(255,0,0), 1, 8, 0);
		cvLine(srcImage, pts[0], pts[2], CV_RGB(255,0,0), 1, 8, 0);

		cvLine(dstImage, pts2[0], pts2[1], CV_RGB(255,0,0), 1, 8, 0);
		cvLine(dstImage, pts2[1], pts2[2], CV_RGB(255,0,0), 1, 8, 0);
		cvLine(dstImage, pts2[0], pts2[2], CV_RGB(255,0,0), 1, 8, 0);
		*tris++;	
	}	
	cvSaveImage("d:\\src.bmp", srcImage);
	cvSaveImage("d:\\dst.bmp", dstImage);

	//	
	AffineTransformFloat(srtPt, dstPt, np, ap);
	ht = srcImage->height;
	wd = srcImage->width;	
	minX = 0;
	minY = 0;
	maxX = wd;
	maxY = ht;	
	cx[0] = 0;  cy[0] = 0;
	cx[1] = wd; cy[1] = 0;
	cx[2] = 0;  cy[2] = ht;
	cx[3] = wd; cy[3] = ht;	
	//calculate the polygon of src image in the dst image
	for(i=0; i<4; i++)
	{
		tx = AFFINE_X(cx[i], cy[i], ap);
		ty = AFFINE_Y(cx[i], cy[i], ap);		
		minX = min(minX, tx);
		minY = min(minY, ty);
		maxX = max(maxX, tx);
		maxY = max(maxY, ty);
	}
	mosaicImage = cvCreateImage( cvSize(maxX-minX, maxY-minY), IPL_DEPTH_8U, 3);
	

	
	printf("\n Mosaic..... \n");
	for(j=0; j<srcImage->height; j++)
	{
		for(i=0; i<srcImage->width; i++)
		{
			//i = 382;
			//j = 308;

			s=cvGet2D(srcImage, j, i); // get the (i,j) pixel value
			
			//find the triangle
			bool bIsfind = false;
			MyPointF srctps[3];
			MyPointF dsttps[3];
			int ti[3];
			MyPointF cp;
			int  minLen = 100000;
			int  len;
			tris = pTin->SelectTriangles(&nTriangleNum, 0, 0, 0, 0 );	
			//printf("\n triangle number: %d \n", nTriangleNum);	
			for(m=0; m<nTriangleNum; m++)
			{		
				for( k=0; k<3; k++)
				{	
					px[k] = (*tris)->vertex[k]->x;
					py[k] = (*tris)->vertex[k]->y;
					ti[k] = (*tris)->vertex[k]->attr;
				}
				
				if(pnpoly(3, px, py, i, j))
				{
					for( k=0; k<3; k++)
					{
						srctps[k].x = srtPt[ ti[k] ].x;
						srctps[k].y = srtPt[ ti[k] ].y;
						dsttps[k].x = dstPt[ ti[k] ].x;
						dsttps[k].y = dstPt[ ti[k] ].y;
					}
					bIsfind = true;
					break;
				}
				else
				{
					cp.x = 0;
					cp.y = 0;
					for( k=0; k<3; k++)
					{
						cp.x += px[k];
						cp.y += py[k];
					}
					cp.x /= 3;
					cp.y /= 3;
					len = sqrt( (cp.x-i)*(cp.x-i) + (cp.y-j)*(cp.y-j) );
					if(len<minLen)
					{
						for( k=0; k<3; k++)
						{
							srctps[k].x = srtPt[ ti[k] ].x;
							srctps[k].y = srtPt[ ti[k] ].y;
							dsttps[k].x = dstPt[ ti[k] ].x;
							dsttps[k].y = dstPt[ ti[k] ].y;
						}
						minLen = len;
					}
				}
				*tris++;	
			}	

			//generate affine based on triangle
			AffineTransformFloat(srctps, dsttps, 3, ap);
			tx = AFFINE_X(i, j, ap);
			ty = AFFINE_Y(i, j, ap);

			if( (ty-minY)>=0 && (ty-minY)<(maxY-minY)
				&& (tx-minX)>=0 && (tx-minX)<(maxX-minX) )

			cvSet2D(mosaicImage, ty-minY, tx-minX, s);
		}
		printf("%d \n", j);
	}

	//cvSaveImage("d:\\mosaic.bmp", mosaicImage);
	
	delete pTin;

	return mosaicImage;
}
*/

/*
   function: 根据两张图片的初始匹配的结果来进行特征点加密
 */
/*
void  AddFeaturePoint(IplImage* srcImage, IplImage* dstImage, 
					  MyPointF** srcPt, MyPointF** dstPt, int* np
					  )
{
	int i,j,k;
	CvScalar cs;
	float ap[6];
	MyPoint* srcFeatpts = new MyPoint[MAX_FEAT_NUM];
	
	MyPoint* srcFeatptsMatching = new MyPoint[MAX_FEAT_NUM];
	MyPoint* dstFeatptsMatching = new MyPoint[MAX_FEAT_NUM];

	float*   sim = new float[MAX_FEAT_NUM];
	int sum;
	int ht,wd;
	unsigned char *usrcImage =  NULL;
	unsigned char *udstImage =  NULL;
	MyPoint curPt;
	float fsim;
	int index;

	ht = dstImage->height;
	wd = dstImage->width;

	//convert from iplimage to unsigned char
	usrcImage = (unsigned char*)malloc(srcImage->height*srcImage->width);
	for(j=0; j<srcImage->height; j++)
		for(i=0; i<srcImage->width; i++)
		{
			cs = cvGet2D(srcImage, j, i);
			usrcImage[j*srcImage->width + i] = cs.val[0];
		}
	SaveBmp("d:\\usrcimage.bmp", usrcImage, srcImage->height, srcImage->width);
	udstImage = (unsigned char*)malloc(dstImage->height*dstImage->width);
	
	for(j=0; j<dstImage->height; j++)
		for(i=0; i<dstImage->width; i++)
		{
			cs = cvGet2D(dstImage, j, i);
			udstImage[j*dstImage->width + i] = cs.val[0];
		}
	SaveBmp("d:\\udstimage.bmp", usrcImage, srcImage->height, srcImage->width);
	
	//建立整体的仿射关系
	AffineTransformFloat(*srcPt, *dstPt, *np, ap);

	//从srcimage中提取特征点
	//DetectFeaturePt(featpts, &sum, uimage, zoomImg1->height, zoomImg1->width);	
	DetectFeaturePt_block(srcFeatpts, &sum, usrcImage, srcImage->height, srcImage->width);	

	
	//基于仿射关系，利用相关系数来进行匹配
	index = 0;
	for(i=0; i<sum; i++)
	{
		MyPoint mp;

		curPt.x = AFFINE_X( srcFeatpts[i].x, srcFeatpts[i].y, ap );		
		curPt.y = AFFINE_Y( srcFeatpts[i].x, srcFeatpts[i].y, ap );	
		
		if( curPt.x<0 || curPt.x>=dstImage->width
			|| curPt.y<0 || curPt.y>=dstImage->height)
			continue;

		mp = PtBlockMatching_Cov_WithRadius(srcFeatpts[i], usrcImage, curPt, udstImage, 
			ht, wd, &fsim, 20 );
		if(fsim>0.8)
		{
			srcFeatptsMatching[index] = srcFeatpts[i];
			dstFeatptsMatching[index] = mp;
			sim[index] = fsim;
			index++;
		}
	}
	
	for(i=0; i<index; i++)
	{
		CvPoint p;
		p.x = dstFeatptsMatching[i].x;
		p.y = dstFeatptsMatching[i].y;
		cvDrawCircle(dstImage, p, 2, CV_RGB(255,0,0));

		p.x = srcFeatptsMatching[i].x;
		p.y = srcFeatptsMatching[i].y;
		cvDrawCircle(srcImage, p, 2, CV_RGB(255,0,0));
	}
	cvSaveImage("d:\\harris-dst.bmp", dstImage);
	cvSaveImage("d:\\harris-src.bmp", srcImage);

	//将加密点与sift点合并
	for(i=0; i<*np; i++)
	{
		dstFeatptsMatching[i+index].x = (*dstPt)[i].x;
		dstFeatptsMatching[i+index].y = (*dstPt)[i].y;		
		srcFeatptsMatching[i+index].x = (*srcPt)[i].x;
		srcFeatptsMatching[i+index].y = (*srcPt)[i].y;
	}
	free(*dstPt);
	free(*srcPt);
	*np = *np + index;
	*dstPt = (MyPointF*)malloc((*np)*sizeof(MyPointF));
	*srcPt = (MyPointF*)malloc((*np)*sizeof(MyPointF));
	for(i=0; i<*np; i++)
	{
		(*dstPt)[i].x = dstFeatptsMatching[i].x;
		(*dstPt)[i].y = dstFeatptsMatching[i].y;		
		(*srcPt)[i].x = srcFeatptsMatching[i].x;
		(*srcPt)[i].y = srcFeatptsMatching[i].y;
	}

	delete[] srcFeatpts;
	delete[] dstFeatptsMatching;
	delete[] srcFeatptsMatching;
	delete[] sim;
	free(usrcImage);
	free(udstImage);
}


IplImage* MosaicUsingHomography(IplImage* srcImage, IplImage* dstImage, double* H, int* nminx, int* nminy)
{
	int i,j,k;
	IplImage* mosaicImage = NULL;
	CvScalar s;
	int tx,ty;
    int cx[4];
	int cy[4];
	double iH[9];

	int ht,wd;
	int minX, minY, maxX, maxY;
	int minX1, minY1, maxX1, maxY1;

	ht = srcImage->height;
	wd = srcImage->width;

	minX = 0;
	minY = 0;
	maxX = dstImage->width;
	maxY = dstImage->height;

	cx[0] = 0;  cy[0] = 0;
	cx[1] = wd; cy[1] = 0;
	cx[2] = 0;  cy[2] = ht;
	cx[3] = wd; cy[3] = ht;

	//calculate the polygon of src image in the dst image
	for(i=0; i<4; i++)
	{
		tx = HOMOGRAPHY_X(cx[i], cy[i], H);
		ty = HOMOGRAPHY_Y(cx[i], cy[i], H);		

		//tx = (H[0]*cx[i] + H[1]*cy[i] + H[2]) / (H[6]*cx[i] + H[7]*cy[i] + H[8]);
		//ty = (H[3]*cx[i] + H[4]*cy[i] + H[5]) / (H[6]*cx[i] + H[7]*cy[i] + H[8]);

		minX = min(minX, tx);
		minY = min(minY, ty);
		maxX = max(maxX, tx);
		maxY = max(maxY, ty);
		cx[i] = tx;
		cy[i] = ty;
	}

	if(nminx!=NULL && nminy!=NULL)
	{
		*nminx = minX;
		*nminy = minY;
	}

	mosaicImage = cvCreateImage( cvSize(maxX-minX, maxY-minY), IPL_DEPTH_8U, 3);

	for(j=0; j<dstImage->height; j++)
		for(i=0; i<dstImage->width; i++)
		{
			s=cvGet2D(dstImage, j, i);
			cvSet2D(mosaicImage, j-minY, i-minX, s);
		}    



	//间接投影的方法
	memcpy(iH, H, 9*sizeof(double));
	invers_matrix(iH, 3);
	minX1 = wd;
	maxX1 = 0;
	minY1 = ht;
	maxY1 = 0;
	for(i=0; i<4; i++)
	{
		if(minX1>cx[i])		minX1 = cx[i];
		if(maxX1<cx[i])		maxX1 = cx[i];
		if(minY1>cy[i])		minY1 = cy[i];
		if(maxY1<cy[i])		maxY1 = cy[i];
	}

	for(j=minY1; j<maxY1; j++)
		for(i=minX1; i<maxX1; i++)
		{
			int ix,iy;
			ix = HOMOGRAPHY_X(i,j,iH);
			iy = HOMOGRAPHY_Y(i,j,iH);
			if(ix>0 && ix<wd && iy>0 && iy<ht)
			{
				s=cvGet2D(srcImage, iy, ix); // get the (i,j) pixel value				
				cvSet2D(mosaicImage, j-minY, i-minX, s);
			}
		}

	return mosaicImage;	
}
*/


#define MAX_IMAGE_SEQ 2000 

//bundler ground to image model
void GrdToImg1(double* R, double* T, double f, double k1, double k2, 
			   double gx, double gy, double gz, double* ix, double* iy, int ht, int wd)
{
	double res[3];
	double P[3];
	double px,py;
	double rp;

	P[0] = gx;
	P[1] = gy;
	P[2] = gz;

	mult(R, P, res, 3, 3, 1);
	res[0] += T[0];
	res[1] += T[1];
	res[2] += T[2];

	px = -res[0]/res[2];
	py = -res[1]/res[2];			

	rp = 1 + k1*(px*px+py*py) + k2*(px*px+py*py)*(px*px+py*py);

	(*ix) = px*f*rp + wd/2;
	(*iy) = py*f*rp + ht/2;
}


//x=RX+T, when the Z is known
void ImgtoGrd1(double* R, double* T, double f, double ix, double iy, double gz, double* gx, double* gy)
{	
	double iR[9];
	double ip[3];
	double iT[3];
	double res[3];

	memcpy(iR, R, sizeof(double)*9);	
	invers_matrix(iR, 3);

	ip[0] = ix;
	ip[1] = iy;
	ip[2] = -f;   
	mult(iR, ip, res, 3, 3, 1);
	mult(iR, T, iT, 3, 3, 1);

	*gx = -res[0]/res[2]*(iT[2]+gz) - iT[0];
	*gy = -res[1]/res[2]*(iT[2]+gz) - iT[1];	
}


void FitPlane1(double* px, double* py, double* pz, int np, double* R)
{	
	int i,j,k;

	// 拟合平面，旋转到与XZ平面平行
	double *errorinfo = new double[np];
	double *errorpoints = new double[np];
	double *tmpotx = new double[np];
	double *tmpoty = new double[np];
	double *tmpotz = new double[np];
	double para[100];
	int paranum = 100;
	int infonum = 0;
	FILE* hp=NULL;
	char  tmppath[256] = "d:\\R.txt";

	FitObject(2,px,py,pz,np,
		para,paranum,
		errorinfo,errorpoints,
		tmpotx,tmpoty,tmpotz,infonum);

	//写文件	
	hp = fopen(tmppath,"wt");
	for (i=0;i<4;i++)
	{
		fprintf(hp, "%d %lf %lf %lf\n", i, tmpotx[i],tmpoty[i],tmpotz[i]);
	}
	fprintf(hp, "%lf %lf %lf %lf\n",para[0],para[1],para[2],para[3]);

	int pairarray[100];
	pairarray[0] = 31;
	int objectinfo[3] = {1,0,0};
	Normalize(para,3);
	double tx[1] = {para[0]};
	double ty[1] = {para[1]};
	double tz[1] = {para[2]};
	ComputeParameter(tx, ty, tz, pairarray, objectinfo, 1, para, paranum);
	fprintf(hp, "%lf %lf %lf\n  %lf %lf %lf\n %lf %lf %lf\n", 
		para[0],para[1],para[2],para[3],para[4],para[5],para[6],para[7],para[8]);
	fclose(hp);

	//rotation matrix from original data to horizontal plane
	for(i=0; i<9; i++)
		R[i] = para[i];	

	delete []errorinfo;
	delete []errorpoints;
	delete []tmpotx;
	delete []tmpoty;
	delete []tmpotz;
}


int  BundlerMosaicOnebyOne(char* listfile, char* bundleFile, int ht, int wd, const char* savePath)
{
	int  numCamera,numPt;
	char strLine[256];
	FILE* fp = NULL;
	int   i,j,k;
	stPOS* pPosSeq=NULL;	
	double gP[3];
	int    rgb[3];
	int    nPt;
	int    nImageIndex,nPtIndex;
	double ix,iy;
	double p[3];
	double res[3];
	double px,py;
	double k1,k2;
	double f;
	double rp;
	double rx,ry;
	double gz = -0.01;
	double bx[4],by[4];
	double gx,gy;
	double minx,maxx,miny,maxy;
	double resolution;
	int    oht,owd;
	IplImage* mosaicImg;
	IplImage* pImg;
	CvScalar  color;
	char      filename[256];	
	double    focus;
	int       it;
	//IplImage** pImgSeq = NULL;
	double* ppx;
	double* ppy;
	double* ppz;
	double  planeR[9];

	bx[0] = -wd/2;	by[0] = -ht/2;
	bx[1] = wd/2;	by[1] = -ht/2;
	bx[2] = -wd/2;	by[2] = ht/2;
	bx[3] = wd/2;	by[3] = ht/2;

	pPosSeq = (stPOS*)malloc(MAX_IMAGE_SEQ*sizeof(stPOS));

	fp =fopen(bundleFile, "r");
	fgets(strLine,256, fp);
	printf("%s\n", strLine);
	fscanf(fp, "%d %d", &numCamera, &numPt);
	printf("%d %d\n", numCamera, numPt);
	for(i=0; i<numCamera; i++)
	{
		fscanf(fp, "%lf %lf %lf", &(pPosSeq[i].f), &(pPosSeq[i].k1), &(pPosSeq[i].k2));
		for(j=0; j<9; j++)	
			fscanf(fp, "%lf", &(pPosSeq[i].R[j]));
		for(j=0; j<3; j++)	
			fscanf(fp, "%lf", &(pPosSeq[i].T[j]));
	}

	//read 3d point
	ppx = (double*)malloc(numPt*sizeof(double));
	ppy = (double*)malloc(numPt*sizeof(double));
	ppz = (double*)malloc(numPt*sizeof(double));
	double *dem=new double[numCamera];
	int *numofpts=new int[numCamera];
	for (i=0;i<numCamera;i++)
	{
		dem[i]=0;
		numofpts[i]=0;
	}
	for(i=0; i<numPt; i++)
	{
		fscanf(fp, "%lf %lf %lf", &gP[0],&gP[1],&gP[2]);
		fscanf(fp, "%d %d %d", &rgb[0],&rgb[1],&rgb[2]);
		fscanf(fp, "%d", &nPt);
		for(k=0; k<nPt; k++)
		{
			fscanf(fp,"%d %d %lf %lf", &nImageIndex, &nPtIndex, &ix, &iy);
			dem[nImageIndex]+=gP[2];
			numofpts[nImageIndex]++;
		}
		ppx[i] = gP[0];
		ppy[i] = gP[1];
		ppz[i] = gP[2];
	}
	fclose(fp);

	FitPlane1(ppx,ppy,ppz,numPt,planeR);
	invers_matrix(planeR,3);
	for (i=0;i<numCamera;i++)
	{
		dem[i]=dem[i]/(numofpts[i]+1);
	}

	double mdem = 0;
	for (i=0;i<numCamera;i++)
	{
		mdem += dem[i];
	}
	mdem /= numCamera;

	/*
	//image dodging, added by xiedonghai, 2012.11.21
	for(k=1; k<numCamera; k++)
	{
	HistMatchingColor(pImgSeq[0], pImgSeq[k]);
	}*/

	//calculate the area of each image
	gz = mdem;
	for(i=0; i<numCamera; i++)
	{
		minx = 100000000;
		maxx = -10000000;
		miny = 100000000;
		maxy = -10000000;
		for(j=0; j<4; j++)
		{
			ix = bx[j];
			iy = by[j];
			ImgtoGrd1(pPosSeq[i].R, pPosSeq[i].T, pPosSeq[i].f, ix, iy, gz, &gx, &gy);

			if(minx>gx) minx=gx;
			if(maxx<gx) maxx=gx;
			if(miny>gy) miny=gy;
			if(maxy<gy) maxy=gy;
		}
		pPosSeq[i].l = minx;
		pPosSeq[i].r = maxx;
		pPosSeq[i].t = miny;
		pPosSeq[i].b = maxy;
	}

	//calculate the ground area of all images
	minx = 100000000;
	maxx = -10000000;
	miny = 100000000;
	maxy = -10000000;
	for(i=0; i<numCamera; i++)
	{
		if(minx>pPosSeq[i].l) minx = pPosSeq[i].l;
		if(maxx<pPosSeq[i].r) maxx = pPosSeq[i].r;
		if(miny>pPosSeq[i].t) miny = pPosSeq[i].t;
		if(maxy<pPosSeq[i].b) maxy = pPosSeq[i].b;
	}
	printf("%lf %lf %lf %lf \n", minx,maxx,miny,maxy);

	//calculate the resolution
	resolution = (pPosSeq[0].r - pPosSeq[0].l) / wd * 1.1;	
	oht = (maxy-miny)/resolution;
	owd = (maxx-minx)/resolution;
    
	double size = 3*oht*owd/1024/1024;

	if(size<480)
	{
		mosaicImg = cvCreateImage(cvSize(owd, oht), IPL_DEPTH_8U, 3);

		fp = fopen(listfile, "r");
		for( k=0; k<numCamera; k++)
		{			
			fscanf(fp, "%s %d %lf", filename, &it, &focus);
			IplImage* pImage = cvLoadImage(filename, 1);
			IplImage* pResizeImage = cvCreateImage( cvSize(wd, ht), pImage->depth, pImage->nChannels);
			cvResize(pImage, pResizeImage);

			for( gy=pPosSeq[k].t; gy<pPosSeq[k].b; gy+=resolution )
			{
				for( gx=pPosSeq[k].l; gx<pPosSeq[k].r; gx+=resolution )
				{
					j = (gy-miny) / resolution;
					i = (gx-minx) / resolution;				
					j = max(0, min(oht-1, j));
					i = max(0, min(owd-1, i));

					double tp[3];
					gP[0] = gx; 
					gP[1] = gy; 
					gP[2] = mdem;
					mult(planeR,gP,tp,3,3,1); 
					double gx1,gy1,gz1;
					gx1 = tp[0];
					gy1 = tp[1];
					gz1 = tp[2];	

					GrdToImg1(pPosSeq[k].R,  pPosSeq[k].T,  pPosSeq[k].f, 
						pPosSeq[k].k1, pPosSeq[k].k2, gx1, gy1, gz1, &ix, &iy, ht, wd);

					if(ix>0 && ix<wd && iy>0 && iy<ht)
					{
						color = cvGet2D(pResizeImage, ht-iy-1, ix);
						cvSet2D(mosaicImg, j, i, color);
					}				
				}
			}
			cvReleaseImage(&pImage);
			cvReleaseImage(&pResizeImage);
		}
		fclose(fp);

		if(savePath!=NULL)
			cvSaveImage(savePath, mosaicImg);
		cvReleaseImage(&mosaicImg);   
	}
	else
	{
		return 0;
	}

	free(ppx);
	free(ppy);
	free(ppz);
	free(dem);
	free(numofpts);
	free(pPosSeq);

	printf("Mosaic finished ! \n");

	return 1;
}

//mosaic based on input files
int  BundleMosaic(char* listfile, char* bundleFile, int ht, int wd, const char* savePath)
{
	int  numCamera,numPt;
	char strLine[256];
	FILE* fp = NULL;
	int   i,j,k;
	stPOS* pPosSeq=NULL;	
	double gP[3];
	int    rgb[3];
	int    nPt;
	int    nImageIndex,nPtIndex;
	double ix,iy;
	double p[3];
	double res[3];
	double px,py;
	double k1,k2;
	double f;
	double rp;
	double rx,ry;
	double gz = -0.01;
	//int    ht,wd;
	double bx[4],by[4];
	double gx,gy;
	double minx,maxx,miny,maxy;
	double resolution;
	int    oht,owd;
	IplImage* mosaicImg;
	IplImage* pImg;
	CvScalar  color;
	char      filename[256];	
	double    focus;
	int       it;
	IplImage** pImgSeq = NULL;
	double* ppx;
	double* ppy;
	double* ppz;
	double  planeR[9];

	bx[0] = -wd/2;	by[0] = -ht/2;
	bx[1] = wd/2;	by[1] = -ht/2;
	bx[2] = -wd/2;	by[2] = ht/2;
	bx[3] = wd/2;	by[3] = ht/2;

	pPosSeq = (stPOS*)malloc(MAX_IMAGE_SEQ*sizeof(stPOS));

	fp =fopen(bundleFile, "r");
	fgets(strLine,256, fp);
	printf("%s\n", strLine);
	fscanf(fp, "%d %d", &numCamera, &numPt);
	printf("%d %d\n", numCamera, numPt);
	for(i=0; i<numCamera; i++)
	{
		fscanf(fp, "%lf %lf %lf", &(pPosSeq[i].f), &(pPosSeq[i].k1), &(pPosSeq[i].k2));
		for(j=0; j<9; j++)	
			fscanf(fp, "%lf", &(pPosSeq[i].R[j]));
		for(j=0; j<3; j++)	
			fscanf(fp, "%lf", &(pPosSeq[i].T[j]));
	}

	//read 3d point
	ppx = (double*)malloc(numPt*sizeof(double));
	ppy = (double*)malloc(numPt*sizeof(double));
	ppz = (double*)malloc(numPt*sizeof(double));
	double *dem=new double[numCamera];
	int *numofpts=new int[numCamera];
	for (i=0;i<numCamera;i++)
	{
		dem[i]=0;
		numofpts[i]=0;
	}
	for(i=0; i<numPt; i++)
	{
		fscanf(fp, "%lf %lf %lf", &gP[0],&gP[1],&gP[2]);
		fscanf(fp, "%d %d %d", &rgb[0],&rgb[1],&rgb[2]);
		fscanf(fp, "%d", &nPt);
		for(k=0; k<nPt; k++)
		{
			fscanf(fp,"%d %d %lf %lf", &nImageIndex, &nPtIndex, &ix, &iy);
			dem[nImageIndex]+=gP[2];
			numofpts[nImageIndex]++;
		}
		ppx[i] = gP[0];
		ppy[i] = gP[1];
		ppz[i] = gP[2];
	}
	fclose(fp);

	FitPlane1(ppx,ppy,ppz,numPt,planeR);
	invers_matrix(planeR,3);
	for (i=0;i<numCamera;i++)
	{
		dem[i]=dem[i]/numofpts[i];
	}

	double mdem = 0;
	for (i=0;i<numCamera;i++)
	{
		mdem += dem[i];
	}
	mdem /= numCamera;

	//loading image
	pImgSeq = (IplImage**)malloc(numCamera*sizeof(IplImage*));	
	fp = fopen(listfile, "r");
	for( k=0; k<numCamera; k++)
	{
		fscanf(fp, "%s %d %lf", filename, &it, &focus);
		//pImgSeq[k] = cvLoadImage( filename, 1);		

		IplImage* pImage = cvLoadImage(filename, 1);
		IplImage* pResizeImage = cvCreateImage( cvSize(wd, ht), pImage->depth, pImage->nChannels);
		cvResize(pImage, pResizeImage);
		pImgSeq[k] = pResizeImage;
		cvReleaseImage(&pImage);
	}
	fclose(fp);

	/*
	//image dodging, added by xiedonghai, 2012.11.21
	for(k=1; k<numCamera; k++)
	{
		HistMatchingColor(pImgSeq[0], pImgSeq[k]);
	}*/
	
	//calculate the area of each image
	gz = mdem;
	for(i=0; i<numCamera; i++)
	{
		minx = 100000000;
		maxx = 0;
		miny = 100000000;
		maxy = 0;
		for(j=0; j<4; j++)
		{
			ix = bx[j];
			iy = by[j];
			ImgtoGrd1(pPosSeq[i].R, pPosSeq[i].T, pPosSeq[i].f, ix, iy, gz, &gx, &gy);

 			if(minx>gx) minx=gx;
			if(maxx<gx) maxx=gx;
			if(miny>gy) miny=gy;
			if(maxy<gy) maxy=gy;
		}
		pPosSeq[i].l = minx;
		pPosSeq[i].r = maxx;
		pPosSeq[i].t = miny;
		pPosSeq[i].b = maxy;
	}

	//calculate the ground area of all images
	minx = 100000000;
	maxx = 0;
	miny = 100000000;
	maxy = 0;
	for(i=0; i<numCamera; i++)
	{
		if(minx>pPosSeq[i].l) minx = pPosSeq[i].l;
		if(maxx<pPosSeq[i].r) maxx = pPosSeq[i].r;
		if(miny>pPosSeq[i].t) miny = pPosSeq[i].t;
		if(maxy<pPosSeq[i].b) maxy = pPosSeq[i].b;
	}
	printf("%lf %lf %lf %lf \n", minx,maxx,miny,maxy);

	//calculate the resolution
	//resolution = 0.002;
	resolution = (pPosSeq[0].r - pPosSeq[0].l) / wd;	
	oht = (maxy-miny)/resolution;
	owd = (maxx-minx)/resolution;

	mosaicImg = cvCreateImage(cvSize(owd,oht), IPL_DEPTH_8U, 3);
	
	for(j=0; j<oht; j++)
	{
		for(i=0; i<owd; i++)
		{
			double tp[3];
			gx = gP[0] = (double)(i)*resolution+minx;
			gy = gP[1] = (double)(j)*resolution+miny;
			gz = gP[2] = mdem;
						
			mult(planeR,gP,tp,3,3,1); 
			
			for( k=0; k<numCamera; k++)
			{
				if( gx<pPosSeq[k].l || gx>pPosSeq[k].r || gy<pPosSeq[k].t || gy>pPosSeq[k].b )
					continue;

				gx = tp[0];
				gy = tp[1];
				gz = tp[2];	

				GrdToImg1(pPosSeq[k].R,  pPosSeq[k].T,  pPosSeq[k].f, 
					pPosSeq[k].k1, pPosSeq[k].k2, gx, gy, gz, &ix, &iy, ht, wd);

				if(ix>0 && ix<wd && iy>0 && iy<ht)
				{
					color = cvGet2D(pImgSeq[k], ht-iy-1, ix);
					cvSet2D(mosaicImg, j, i, color);
				}
			}
		}
		printf("%d \n", j);
	}

	if(savePath!=NULL)
		cvSaveImage(savePath, mosaicImg);
	cvReleaseImage(&mosaicImg);    

	for(i=0; i<numCamera; i++)
	{
		cvReleaseImage( &pImgSeq[i] );
	}
	free(ppx);
	free(ppy);
	free(ppz);
	free(pImgSeq);
	free(dem);
	free(numofpts);
    free(pPosSeq);

	printf("Mosaic finished ! \n");
	return 1;
}

int PlaneMosaic(char* listfile, char* bundleFile, double ratio, const char* savePath)
{

	int  numCamera,numPt;
	char strLine[256];
	FILE* fp = NULL;
	int   i,j,k,t;
	stPOS* pPosSeq=NULL;	
	double gP[3];
	int    rgb[3];
	int    nPt;
	int    nImageIndex,nPtIndex;
	double ix,iy;
	double p[3];
	double res[3];
	double px,py;
	double k1,k2;
	double f;
	double rp;
	double rx,ry;
	double gz = -0.01;
	double bx[4],by[4];
	double gx,gy;
	double minx,maxx,miny,maxy;
	double resolution;
	int    oht,owd;
	IplImage* mosaicImg;
	IplImage* pImg;
	CvScalar  color;
	char      filename[256];	
	double    focus;
	int       it;
	//IplImage** pImgSeq = NULL;
	double* ppx;
	double* ppy;
	double* ppz;
	double  planeR[9];


	//reading list file
	char imagefile[256];
	fp = fopen(listfile, "r");
	fscanf(fp, "%s %d %lf ", imagefile, &t, &focus);
	fclose(fp);

	IplImage* pImage = cvLoadImage(imagefile);
	int ht = pImage->height;
	int wd = pImage->width;
    cvReleaseImage(&pImage); 

	bx[0] = -wd/2;	by[0] = -ht/2;
	bx[1] = wd/2;	by[1] = -ht/2;
	bx[2] = -wd/2;	by[2] = ht/2;
	bx[3] = wd/2;	by[3] = ht/2;

	pPosSeq = (stPOS*)malloc(MAX_IMAGE_SEQ*sizeof(stPOS));

	fp =fopen(bundleFile, "r");
	fgets(strLine,256, fp);
	printf("%s\n", strLine);
	fscanf(fp, "%d %d", &numCamera, &numPt);
	printf("%d %d\n", numCamera, numPt);
	for(i=0; i<numCamera; i++)
	{
		fscanf(fp, "%lf %lf %lf", &(pPosSeq[i].f), &(pPosSeq[i].k1), &(pPosSeq[i].k2));
		for(j=0; j<9; j++)	
			fscanf(fp, "%lf", &(pPosSeq[i].R[j]));
		for(j=0; j<3; j++)	
			fscanf(fp, "%lf", &(pPosSeq[i].T[j]));
	}

	//read 3d point
	ppx = (double*)malloc(numPt*sizeof(double));
	ppy = (double*)malloc(numPt*sizeof(double));
	ppz = (double*)malloc(numPt*sizeof(double));
	double *dem=new double[numCamera];
	int *numofpts=new int[numCamera];
	for (i=0;i<numCamera;i++)
	{
		dem[i]=0;
		numofpts[i]=0;
	}
	for(i=0; i<numPt; i++)
	{
		fscanf(fp, "%lf %lf %lf", &gP[0],&gP[1],&gP[2]);
		fscanf(fp, "%d %d %d", &rgb[0],&rgb[1],&rgb[2]);
		fscanf(fp, "%d", &nPt);
		for(k=0; k<nPt; k++)
		{
			fscanf(fp,"%d %d %lf %lf", &nImageIndex, &nPtIndex, &ix, &iy);
			//dem[nImageIndex]+=gP[2];
			numofpts[nImageIndex]++;
		}
		ppx[i] = gP[0];
		ppy[i] = gP[1];
		ppz[i] = gP[2];
	}
	fclose(fp);

	FitPlane1(ppx,ppy,ppz,numPt,planeR);
	
	//transform
	double mdem = 0;
	//fp = fopen("d:\\groud.txt", "w");
	for( i=0; i<numPt; i++)
	{
		double tp[3];
		gP[0] = ppx[i]; 
		gP[1] = ppy[i]; 
		gP[2] = ppz[i];
		mult(planeR,gP,tp,3,3,1);
		ppx[i] = tp[0]; 
		ppy[i] = tp[1]; 
		ppz[i] = tp[2];
		//fprintf(fp, "%lf %lf  \n", gP[2], tp[2]);
		mdem += tp[2];
	}
	//fclose(fp);

	mdem /= numPt;

	//inverse the plane rotation
	invers_matrix(planeR,3);

	/*
	//image dodging, added by xiedonghai, 2012.11.21
	for(k=1; k<numCamera; k++)
	{
	HistMatchingColor(pImgSeq[0], pImgSeq[k]);
	}*/

	//calculate the area of each image
	gz = mdem;
	for(i=0; i<numCamera; i++)
	{
		minx = 100000000;
		maxx = -10000000;
		miny = 100000000;
		maxy = -10000000;
		
		if(pPosSeq[i].f==0)
			continue;

		//transform the free 3d points into ground coordinate
		double tR[9];
		mult(pPosSeq[i].R, planeR, tR, 3, 3, 3);

		for(j=0; j<4; j++)
		{
			ix = bx[j];
			iy = by[j];
			ImgtoGrd1(tR, pPosSeq[i].T, pPosSeq[i].f, ix, iy, gz, &gx, &gy);

			if(minx>gx) minx=gx;
			if(maxx<gx) maxx=gx;
			if(miny>gy) miny=gy;
			if(maxy<gy) maxy=gy;
		}
		pPosSeq[i].l = minx;
		pPosSeq[i].r = maxx;
		pPosSeq[i].t = miny;
		pPosSeq[i].b = maxy;
	}

	//calculate the ground area of all images
	minx = 100000000;
	maxx = -10000000;
	miny = 100000000;
	maxy = -10000000;
	for(i=0; i<numCamera; i++)
	{
		if(pPosSeq[i].f==0)
			continue;
		if(minx>pPosSeq[i].l) minx = pPosSeq[i].l;
		if(maxx<pPosSeq[i].r) maxx = pPosSeq[i].r;
		if(miny>pPosSeq[i].t) miny = pPosSeq[i].t;
		if(maxy<pPosSeq[i].b) maxy = pPosSeq[i].b;
	}
	printf("%lf %lf %lf %lf \n", minx,maxx,miny,maxy);

	//calculate the resolution
	resolution = 0;
	int    nRightCam = 0;
	int    dstWd = wd*ratio;
	for(i=0; i<numCamera; i++)
	{
		if(pPosSeq[i].f==0)
			continue;
		resolution += (pPosSeq[i].r - pPosSeq[i].l) / dstWd;
		nRightCam++;
	}
	resolution = (resolution / nRightCam ) * 1.1;	
	oht = (maxy-miny)/resolution;
	owd = (maxx-minx)/resolution;

	double size = 3*oht*owd/1024/1024;

	//image size is limited to less than 480MB in 32bit OS 
	if(size<480)
	{
		mosaicImg = cvCreateImage(cvSize(owd, oht), IPL_DEPTH_8U, 3);

		fp = fopen(listfile, "r");
		for( k=0; k<numCamera; k++)
		{	
			fscanf(fp, "%s %d %lf", filename, &it, &focus);
			printf("%s \n", filename);

			if(pPosSeq[k].f==0)
				continue;

			IplImage* pImage = cvLoadImage(filename, 1);
			//IplImage* pResizeImage = cvCreateImage( cvSize(wd, ht), pImage->depth, pImage->nChannels);
			//cvResize(pImage, pResizeImage);

			for( gy=pPosSeq[k].t; gy<pPosSeq[k].b; gy+=resolution )
			{
				for( gx=pPosSeq[k].l; gx<pPosSeq[k].r; gx+=resolution )
				{
					j = (gy-miny) / resolution;
					i = (gx-minx) / resolution;				
					j = max(0, min(oht-1, j));
					i = max(0, min(owd-1, i));

					double tp[3];
					gP[0] = gx; 
					gP[1] = gy; 
					gP[2] = mdem;
					mult(planeR,gP,tp,3,3,1); 
					double gx1,gy1,gz1;
					gx1 = tp[0];
					gy1 = tp[1];
					gz1 = tp[2];	

					GrdToImg1(pPosSeq[k].R,  pPosSeq[k].T,  pPosSeq[k].f, 
						pPosSeq[k].k1, pPosSeq[k].k2, gx1, gy1, gz1, &ix, &iy, ht, wd);

					if(ix>0 && ix<wd && iy>0 && iy<ht)
					{
						//color = cvGet2D(pResizeImage, ht-iy-1, ix);
						color = cvGet2D(pImage, ht-iy-1, ix);
						cvSet2D(mosaicImg, j, i, color);
					}				
				}
			}
			cvReleaseImage(&pImage);
			//cvReleaseImage(&pResizeImage);
		}
		fclose(fp);

		if(savePath!=NULL)
			cvSaveImage(savePath, mosaicImg);
		cvReleaseImage(&mosaicImg);   
	}
	else
	{
		return 0;
	}

	free(ppx);
	free(ppy);
	free(ppz);
	free(dem);
	free(numofpts);
	free(pPosSeq);

	printf("Mosaic finished ! \n");

	return 1;
}


//mosaic using absolute orientation
int  BundlerMosaicAR(char* listfile, char* bundleFile, int ht, int wd, const char* savePath)
{
	int  numCamera,numPt;
	char strLine[256];
	FILE* fp = NULL;
	int   i,j,k;
	stPOS* pPosSeq=NULL;	
	double gP[3];
	int    rgb[3];
	int    nPt;
	int    nImageIndex,nPtIndex;
	double ix,iy;
	double p[3];
	double res[3];
	double px,py;
	double k1,k2;
	double f;
	double rp;
	double rx,ry;
	double gz = -0.01;
	double bx[4],by[4];
	double gx,gy;
	double minx,maxx,miny,maxy;
	double resolution;
	int    oht,owd;
	IplImage* mosaicImg;
	IplImage* pImg;
	CvScalar  color;
	char      filename[256];	
	double    focus;
	int       it;
	//IplImage** pImgSeq = NULL;
	double* ppx;
	double* ppy;
	double* ppz;
	double  planeR[9];

	bx[0] = -wd/2;	by[0] = -ht/2;
	bx[1] = wd/2;	by[1] = -ht/2;
	bx[2] = -wd/2;	by[2] = ht/2;
	bx[3] = wd/2;	by[3] = ht/2;

	pPosSeq = (stPOS*)malloc(MAX_IMAGE_SEQ*sizeof(stPOS));

	fp =fopen(bundleFile, "r");
	fgets(strLine,256, fp);
	printf("%s\n", strLine);
	fscanf(fp, "%d %d", &numCamera, &numPt);
	printf("%d %d\n", numCamera, numPt);
	for(i=0; i<numCamera; i++)
	{
		fscanf(fp, "%lf %lf %lf", &(pPosSeq[i].f), &(pPosSeq[i].k1), &(pPosSeq[i].k2));
		for(j=0; j<9; j++)	
			fscanf(fp, "%lf", &(pPosSeq[i].R[j]));
		for(j=0; j<3; j++)	
			fscanf(fp, "%lf", &(pPosSeq[i].T[j]));
	}

	//read 3d point
	ppx = (double*)malloc(numPt*sizeof(double));
	ppy = (double*)malloc(numPt*sizeof(double));
	ppz = (double*)malloc(numPt*sizeof(double));
	double *dem=new double[numCamera];
	int *numofpts=new int[numCamera];
	for (i=0;i<numCamera;i++)
	{
		dem[i]=0;
		numofpts[i]=0;
	}
	for(i=0; i<numPt; i++)
	{
		fscanf(fp, "%lf %lf %lf", &gP[0],&gP[1],&gP[2]);
		fscanf(fp, "%d %d %d", &rgb[0],&rgb[1],&rgb[2]);
		fscanf(fp, "%d", &nPt);
		for(k=0; k<nPt; k++)
		{
			fscanf(fp,"%d %d %lf %lf", &nImageIndex, &nPtIndex, &ix, &iy);
			//dem[nImageIndex]+=gP[2];
			numofpts[nImageIndex]++;
		}
		ppx[i] = gP[0];
		ppy[i] = gP[1];
		ppz[i] = gP[2];
	}
	fclose(fp);

	FitPlane1(ppx,ppy,ppz,numPt,planeR);
	
	//transform
	double mdem = 0;
	//fp = fopen("d:\\groud.txt", "w");
	for( i=0; i<numPt; i++)
	{
		double tp[3];
		gP[0] = ppx[i]; 
		gP[1] = ppy[i]; 
		gP[2] = ppz[i];
		mult(planeR,gP,tp,3,3,1);
		ppx[i] = tp[0]; 
		ppy[i] = tp[1]; 
		ppz[i] = tp[2];
		//fprintf(fp, "%lf %lf  \n", gP[2], tp[2]);
		mdem += tp[2];
	}
	//fclose(fp);

	mdem /= numPt;

	//inverse the plane rotation
	invers_matrix(planeR,3);

	/*
	//image dodging, added by xiedonghai, 2012.11.21
	for(k=1; k<numCamera; k++)
	{
	HistMatchingColor(pImgSeq[0], pImgSeq[k]);
	}*/

	//calculate the area of each image
	gz = mdem;
	for(i=0; i<numCamera; i++)
	{
		minx = 100000000;
		maxx = -10000000;
		miny = 100000000;
		maxy = -10000000;
		
		if(pPosSeq[i].f==0)
			continue;

		//transform the free 3d points into ground coordinate
		double tR[9];
		mult(pPosSeq[i].R, planeR, tR, 3, 3, 3);

		for(j=0; j<4; j++)
		{
			ix = bx[j];
			iy = by[j];
			ImgtoGrd1(tR, pPosSeq[i].T, pPosSeq[i].f, ix, iy, gz, &gx, &gy);

			if(minx>gx) minx=gx;
			if(maxx<gx) maxx=gx;
			if(miny>gy) miny=gy;
			if(maxy<gy) maxy=gy;
		}
		pPosSeq[i].l = minx;
		pPosSeq[i].r = maxx;
		pPosSeq[i].t = miny;
		pPosSeq[i].b = maxy;
	}

	//calculate the ground area of all images
	minx = 100000000;
	maxx = -10000000;
	miny = 100000000;
	maxy = -10000000;
	for(i=0; i<numCamera; i++)
	{
		if(pPosSeq[i].f==0)
			continue;
		if(minx>pPosSeq[i].l) minx = pPosSeq[i].l;
		if(maxx<pPosSeq[i].r) maxx = pPosSeq[i].r;
		if(miny>pPosSeq[i].t) miny = pPosSeq[i].t;
		if(maxy<pPosSeq[i].b) maxy = pPosSeq[i].b;
	}
	printf("%lf %lf %lf %lf \n", minx,maxx,miny,maxy);

	//calculate the resolution
	resolution = 0;
	int    nRightCam = 0;
	for(i=0; i<numCamera; i++)
	{
		if(pPosSeq[i].f==0)
			continue;
		resolution += (pPosSeq[i].r - pPosSeq[i].l) / wd;
		nRightCam++;
	}
	resolution = (resolution / nRightCam ) * 1.1;	
	oht = (maxy-miny)/resolution;
	owd = (maxx-minx)/resolution;

	double size = 3*oht*owd/1024/1024;

	//image size is limited to less than 480MB in 32bit OS 
	if(size<480)
	{
		mosaicImg = cvCreateImage(cvSize(owd, oht), IPL_DEPTH_8U, 3);

		fp = fopen(listfile, "r");
		for( k=0; k<numCamera; k++)
		{	
			fscanf(fp, "%s %d %lf", filename, &it, &focus);

			if(pPosSeq[k].f==0)
				continue;

			IplImage* pImage = cvLoadImage(filename, 1);
			IplImage* pResizeImage = cvCreateImage( cvSize(wd, ht), pImage->depth, pImage->nChannels);
			cvResize(pImage, pResizeImage);

			for( gy=pPosSeq[k].t; gy<pPosSeq[k].b; gy+=resolution )
			{
				for( gx=pPosSeq[k].l; gx<pPosSeq[k].r; gx+=resolution )
				{
					j = (gy-miny) / resolution;
					i = (gx-minx) / resolution;				
					j = max(0, min(oht-1, j));
					i = max(0, min(owd-1, i));

					double tp[3];
					gP[0] = gx; 
					gP[1] = gy; 
					gP[2] = mdem;
					mult(planeR,gP,tp,3,3,1); 
					double gx1,gy1,gz1;
					gx1 = tp[0];
					gy1 = tp[1];
					gz1 = tp[2];	

					GrdToImg1(pPosSeq[k].R,  pPosSeq[k].T,  pPosSeq[k].f, 
						pPosSeq[k].k1, pPosSeq[k].k2, gx1, gy1, gz1, &ix, &iy, ht, wd);

					if(ix>0 && ix<wd && iy>0 && iy<ht)
					{
						color = cvGet2D(pResizeImage, ht-iy-1, ix);
						cvSet2D(mosaicImg, j, i, color);
					}				
				}
			}
			cvReleaseImage(&pImage);
			cvReleaseImage(&pResizeImage);
		}
		fclose(fp);

		if(savePath!=NULL)
			cvSaveImage(savePath, mosaicImg);
		cvReleaseImage(&mosaicImg);   
	}
	else
	{
		return 0;
	}

	free(ppx);
	free(ppy);
	free(ppz);
	free(dem);
	free(numofpts);
	free(pPosSeq);

	printf("Mosaic finished ! \n");

	return 1;
}
