
//corelib
#include "commondata.h"

//
#include "ellipse.h"


#include <vector>
using namespace std;





void DrawEllipse(stEllipse ellipsePara, IplImage* pImage)
{

	int ht = pImage->height;
	int wd = pImage->width;
	int swd = pImage->widthStep;

	int x0 = ellipsePara.x0;
	int y0 = ellipsePara.y0;
	int a = ellipsePara.a;
	int b = ellipsePara.b;

	/*
	int il = min(wd-1, max(0, x0-a));
	int ir = min(wd-1, max(0, x0+a));
	int it = min(ht-1, max(0, y0-b));
	int ib = min(ht-1, max(0, y0+b));
	*/
	
	CvPoint  pts[360];

	for(int i=0; i<360; i++)
	{
		double cosBeita  = cos( -ellipsePara.angle*DPI );
		double sinBeita  = sin( -ellipsePara.angle*DPI );

		double x = a*cos(i*DPI);
		double y = b*sin(i*DPI);

		double nx =  x;//-x0;
		double ny =  y;//-y0;

		double rx = nx*cosBeita + ny*sinBeita;
		double ry = -nx*sinBeita + ny*cosBeita;
         
		rx += x0;
		ry += y0;
        

		pts[i].x = rx;
		pts[i].y = ry;		
	}

	for(int i=0; i<359; i++)
	{
		CvPoint p1,p2;
		p1.x = pts[i].x;
		p1.y = pts[i].y;
		p2.x = pts[i+1].x;
		p2.y = pts[i+1].y;		

		if( p1.x>0 && p1.x<wd && p1.y>0 && p1.y<ht
			&& p2.x>0 && p2.x<wd && p2.y>0 && p2.y<ht)
		 {
		    cvDrawLine(pImage, p1, p2, CV_RGB(255,0,0));
		 }
	}
}

double CalcLongestDis(vector<POINT2> pts)
{
	double dis = 0;

	for(int j=0; j<pts.size(); j++)
		for(int i=j+1; i<pts.size(); i++)
		{
			double d = sqrt( (pts[j].x-pts[i].x)*(pts[j].x-pts[i].x) + (pts[j].y-pts[i].y)*(pts[j].y-pts[i].y)  );
			if(d>dis)
				dis = d;
		}
		return dis;
}

int RandomHoughEllipse(vector<POINT2> edgePts, vector<stEllipse>& vecEllipse, double minDistance, int minCount)
{
	//random select edge pixels
	CvRNG rng = cvRNG(-1);
	int count = edgePts.size();
	//printf("count: %d \n", count);

	if(count<18)
		return 0;

	//double longestDis = CalcLongestDis(edgePts);
	//if( longestDis<18 )
	//	return 0;

	for(int it=0; it<count; it++ )
	{		
		//6 points random selected
		vector<POINT2> randomPts;
		int nRandomNum = 0;
		bool bIsRandomOK = true;
		for(int i=0; i<6; i++)
		{
			nRandomNum++;
			if( nRandomNum>(4*count) )
			{
				bIsRandomOK = false;
				break;
			}

			int idx = cvRandInt(&rng) % count;
			POINT2 pt = edgePts[idx];
			
			int minLen = 10000;
			for(int k=0; k<randomPts.size(); k++)
			{
				POINT2 tp = randomPts[k];
				int len = sqrt( (pt.x-tp.x)*(pt.x-tp.x) + (pt.y-tp.y)*(pt.y-tp.y) );
				if(len<minLen)
					minLen = len;
			}

			if(minLen<2 )
			{
				i--;
				continue;
			}

			randomPts.push_back(pt);
		}

		if(!bIsRandomOK)
			break;

		//ellipse fitting
		stEllipse ellipsePara;
		double err = 0;
		InvertEllipse(randomPts, ellipsePara, err);

		if( err>(minDistance) )
			continue;

		//finding the edge points close to the ellipse
		vector<int> closeIdx;
		for(int i=0; i<edgePts.size(); i++)
		{
			POINT2 pt = edgePts[i];
			double err = Pt2Ellipse(pt.x, pt.y, ellipsePara); //EllipseFunc(pt.x, pt.y, ellipsePara);
			if( err<minDistance)
			{
				closeIdx.push_back(i);
			}
		}

		if( closeIdx.size()>minCount )// && fabs(ellipsePara.angle)<15)
		{			
			vecEllipse.push_back(ellipsePara);		

			//////  remove the found points //////////
			vector<int> label;
			label.resize(edgePts.size(), 1);
			for(int i=0; i<closeIdx.size(); i++)
			{				
				label[ closeIdx[i] ] = 0;
			}
			vector<POINT2> remainedPts;
			for(int i=0; i<edgePts.size(); i++)
			{
				if( label[i] )
					remainedPts.push_back( edgePts[i] );
			}
			edgePts = remainedPts;
			count = edgePts.size();
			////////////////////////////////////////////

			if(count<12)
				break;
		}	

		//cvSaveImage("d:\\edge.jpg", pInput);
		//cvSaveImage("d:\\ellipse.jpg", pShow);
	}

	return 1;
}

int DetectEllipseUsingContour(IplImage* pInput, vector<stEllipse>& ellipseArray)
{

	//1. finding contours
	CvMemStorage* stor = cvCreateMemStorage(0);
	CvSeq* cont = cvCreateSeq(CV_SEQ_ELTYPE_POINT, sizeof(CvSeq), sizeof(CvPoint) , stor);
	int Nc = cvFindContours(pInput, stor, &cont);

	/*
	//draw contours
	IplImage* pDisp = cvCreateImage( cvGetSize(pInput), 8, 3 );
	memset(pDisp->imageData, 0, pDisp->height*pDisp->widthStep);
	int n=0; 
	printf("Total Contours Detected: %d \n", Nc ); 
	for( CvSeq* c=cont; c!=NULL; c=c->h_next ) 
	{ 
		int r = (double)( rand() ) / (double)(RAND_MAX) * 255;
		int g = (double)( rand() ) / (double)(RAND_MAX) * 255;
		int b = (double)( rand() ) / (double)(RAND_MAX) * 255;
		printf("%d %d %d \n", r,g,b);
		cvDrawContours( pDisp, c, CV_RGB(r,g,b), CV_RGB(r,g,b),0,1,8);
	}
	cvSaveImage("d:\\contours.jpg", pDisp);
	*/

	//2. process each contour
	CvPoint* PointArray;
	CvPoint2D32f* PointArray2D32f;
	for(; cont; cont = cont->h_next)
	{   
		printf(".");

		int i; // Indicator of cycle.
		int count = cont->total; // This is number point in contour
		CvPoint center;
		CvSize size;

		// Number point must be more than or equal to 6 (for cvFitEllipse_32f).        
		if( count < 6 )
			continue;

		// Alloc memory for contour point set.    
		PointArray = (CvPoint*)malloc( count*sizeof(CvPoint) );
		// Get contour point set.
		cvCvtSeqToArray(cont, PointArray, CV_WHOLE_SEQ);

		// Convert CvPoint set to CvBox2D32f set.
		vector<POINT2> pts;
		for(i=0; i<count; i++)
		{
			POINT2 tp;
			tp.x = (float)PointArray[i].x;
			tp.y = (float)PointArray[i].y;
			pts.push_back(tp);	
		}
		free(PointArray);
		
		vector<stEllipse> vecEllipse;
		RandomHoughEllipse(pts, vecEllipse, 0.6, 30);

		for(int i=0; i<vecEllipse.size(); i++)
		{
			ellipseArray.push_back( vecEllipse[i] );
		}
	}
	
	return 1;
}

//detect ellipse using all edge points of the image
int DetectEllipse(IplImage* pInput, vector<stEllipse>& ellipseArray)
{	
	double minDistance = 0.6;  //the distance from the point to ellipse 
	double minCount = 64;    //the minimal pixels of the ellipse
    
	//
	vector<POINT2> edgePts;

	int ht = pInput->height;
	int wd = pInput->width;
	int swd = pInput->widthStep;

	//collect the edge pixels
	for(int j=0; j<ht; j++)
		for(int i=0; i<wd; i++)
		{
			unsigned char value = (unsigned char)(pInput->imageData[j*swd+i]);
			if(value>0)
			{
				POINT2 pt;
				pt.x = i;
				pt.y = j;
				edgePts.push_back(pt);
			}
		}

	vector<int> label;
	label.resize(edgePts.size());
	for(int i=0; i<label.size(); i++)
		label[i] = 1;

	//random select edge pixels
	CvRNG rng = cvRNG(-1);
	int count = edgePts.size();
	int nMarkPt = 0;

	printf("count: %d \n", count);

	for(int it=0; it<count; it++ )
	{		
		printf("%d \n", it);

		//6 points random selected
		vector<POINT2> randomPts;
		for(int i=0; i<6; i++)
		{
			 int idx = cvRandInt(&rng) % count;
			 POINT2 pt = edgePts[idx];

			 /*
			 int value = (unsigned char)(pInput->imageData[ (int)(pt.y)*swd+(int)(pt.x)]);
			 if(value==0)
			 {
				 i--;
				 continue;
			 }*/
             
			 
			 int bIsRemove = 1-label[idx];
			 if(bIsRemove>0)
			 {   
				 i--;
				 continue;
			 }
			
			 int minLen = 1000000;
			 for(int k=0; k<randomPts.size(); k++)
			 {
				POINT2 tp = randomPts[k];
				int len = sqrt( (pt.x-tp.x)*(pt.x-tp.x) + (pt.y-tp.y)*(pt.y-tp.y) );
				if(len<minLen)
					minLen = len;
			 }

			 if(minLen<4)
			 {
				 i--;
				 continue;
			 }
			 
			 randomPts.push_back(pt);
		}

		//ellipse fitting
		stEllipse ellipsePara;
		double err = 0;
		InvertEllipse(randomPts, ellipsePara, err);
        
		if(err>minDistance)
			continue;

		if(ellipsePara.a>120)
			continue;

		//finding the edge points close to the ellipse
		//vector<POINT2> closePts;
		vector<int> closeIdx;
		for(int i=0; i<edgePts.size(); i++)
		{
			POINT2 pt = edgePts[i];
		 	double err = Pt2Ellipse(pt.x, pt.y, ellipsePara); //EllipseFunc(pt.x, pt.y, ellipsePara);
			if( err<minDistance)
			{
				//closePts.push_back(pt);
				closeIdx.push_back(i);
			}
		}

		if( closeIdx.size()>minCount )// && fabs(ellipsePara.angle)<15)
		{			
			ellipseArray.push_back(ellipsePara);
			
			/*
			for(int k=0; k<closeIdx.size(); k++)
			{
				label[closeIdx[k]]=0;
			}
			*/

			/*
			for(int k=0; k<closePts.size(); k++)
			{
				int x = closePts[k].x;
				int y = closePts[k].y;				
				//pInput->imageData[y*swd+x] = 0;
			}
			nMarkPt += closePts.size();
			*/
		}

		int nRemains = 0;
		for(int i=0; i<label.size(); i++)
		{
			nRemains += label[i];
		}

		if(nRemains<16)
			break;

		//cvSaveImage("d:\\edge.jpg", pInput);
		//cvSaveImage("d:\\ellipse.jpg", pShow);
	}

	return ellipseArray.size();
}

//calculate the ellipse parameters using more than six points
int InvertEllipse(vector<POINT2> pts, stEllipse& ellipsePara, double& err)
{
	CvBox2D32f* box;
	//CvPoint* PointArray;
	CvPoint2D32f* PointArray2D32f;

	int count = pts.size();

	if(count<6)
		return 0;

	// Alloc memory for contour point set.    
	//PointArray = (CvPoint*)malloc( count*sizeof(CvPoint) );
	PointArray2D32f= (CvPoint2D32f*)malloc( count*sizeof(CvPoint2D32f) );

	// Alloc memory for ellipse data.
	box = (CvBox2D32f*)malloc(sizeof(CvBox2D32f));

	// Get contour point set.
	//cvCvtSeqToArray(cont, PointArray, CV_WHOLE_SEQ);

	// Convert CvPoint set to CvBox2D32f set.
	for(int i=0; i<count; i++)
	{
		PointArray2D32f[i].x = pts[i].x; //(float)PointArray[i].x;
		PointArray2D32f[i].y = pts[i].y; //(float)PointArray[i].y;
	}

	// Fits ellipse to current contour.
	cvFitEllipse(PointArray2D32f, count, box);

	ellipsePara.a  = box->size.height*0.5;
	ellipsePara.b  = box->size.width*0.5;
	ellipsePara.x0	  = box->center.x;
	ellipsePara.y0    = box->center.y;
	ellipsePara.angle = box->angle-90;
	//ellipsePara.angle = 90-box->angle;


	//calculate the error
	//double err = 0;
	for(int i=0; i<count; i++)
	{
		err += Pt2Ellipse(pts[i].x, pts[i].y, ellipsePara);		
	}
	err /= double(count); 

	free(PointArray2D32f);
	free(box);
}

double EllipseFunc(double x, double y, stEllipse ellipsePara)
{
	double res = 0;

	double a = ellipsePara.a;
	double b = ellipsePara.b;

	double cosBeita = cos( ellipsePara.angle*DPI );
	double sinBeita = sin( ellipsePara.angle*DPI );

	double nx =  x-ellipsePara.x0;
	double ny =  y-ellipsePara.y0;
    
	double rx = nx*cosBeita + ny*sinBeita;
	double ry = -nx*sinBeita + ny*cosBeita;

	//double nx = ellipsePara.x0 + (x*cosBeita - y*sinBeita);
	//double ny = ellipsePara.y0 - (x*sinBeita + y*cosBeita);
	
	res = (rx*rx)/(a*a) + (ry*ry)/(b*b) - 1;

	return fabs(res);
}

//the distance from a point to the ellipse
double Pt2Ellipse(double x, double y, stEllipse ellipsePara)
{
	double a = ellipsePara.a;
	double b = ellipsePara.b;

	double cosBeita = cos( ellipsePara.angle*DPI );
	double sinBeita = sin( ellipsePara.angle*DPI );

	double nx =  x-ellipsePara.x0;
	double ny =  y-ellipsePara.y0;

	double rx = nx*cosBeita + ny*sinBeita;
	double ry = -nx*sinBeita + ny*cosBeita;

	double iy = sqrt( (a*a*b*b) / ( a*a+(b*b*rx*rx)/(ry*ry) ) );

	if( ry<0 ) iy = -iy;

	double ix = iy*rx/ry;

	double dis = sqrt( (rx-ix)*(rx-ix) +  (ry-iy)*(ry-iy) );

	return dis;
}




CEllipseDetectOneByOne::CEllipseDetectOneByOne()
{

}

CEllipseDetectOneByOne::~CEllipseDetectOneByOne()
{

}

int CEllipseDetectOneByOne::Detect(IplImage* pInput, vector<stEllipse>& ellipseArray, double fitErr, double minCount)
{
	//1. finding contours
	CvMemStorage* stor = cvCreateMemStorage(0);
	CvSeq* cont = cvCreateSeq(CV_SEQ_ELTYPE_POINT, sizeof(CvSeq), sizeof(CvPoint) , stor);
	int Nc = cvFindContours(pInput, stor, &cont);

	
	//draw contours
	IplImage* pDisp = cvCreateImage( cvGetSize(pInput), 8, 3 );
	memset(pDisp->imageData, 0, pDisp->height*pDisp->widthStep);
	int n=0; 
	printf("Total Contours Detected: %d \n", Nc ); 
	int ic = 1;
	for( CvSeq* c=cont; c!=NULL; c=c->h_next ) 
	{ 
	int r = (double)( rand() ) / (double)(RAND_MAX) * 255;
	int g = (double)( rand() ) / (double)(RAND_MAX) * 255;
	int b = (double)( rand() ) / (double)(RAND_MAX) * 255;
	//int r = ic*20;
	//int g = 0;
	//int b = 0;
	printf("%d %d %d \n", r,g,b);
	cvDrawContours( pDisp, c, CV_RGB(r,g,b), CV_RGB(r,g,b),0,1,8);
	}
	cvSaveImage("d:\\contours.jpg", pDisp);
	

	//2. process each contour
	CvPoint* PointArray;
	CvPoint2D32f* PointArray2D32f;
	for(; cont; cont = cont->h_next)
	{   
		printf(".");

		int i; // Indicator of cycle.
		int count = cont->total; // This is number point in contour
		CvPoint center;
		CvSize size;

		// Number point must be more than or equal to 6 (for cvFitEllipse_32f).        
		if( count < 6 )
			continue;

		// Alloc memory for contour point set.    
		PointArray = (CvPoint*)malloc( count*sizeof(CvPoint) );
		// Get contour point set.
		cvCvtSeqToArray(cont, PointArray, CV_WHOLE_SEQ);
		vector<POINT2> pts;
		for(i=0; i<count; i++)
		{
			POINT2 tp;
			tp.x = (float)PointArray[i].x;
			tp.y = (float)PointArray[i].y;
			pts.push_back(tp);	
		}
		free(PointArray);

		vector<stEllipse> vecEllipse;
		RandomHoughEllipse(pts, vecEllipse, fitErr, minCount);

		for(int i=0; i<vecEllipse.size(); i++)
		{
			ellipseArray.push_back( vecEllipse[i] );
		}
	}

	return 1;
}


CEllipseDetectFromAll::CEllipseDetectFromAll()
{

}

CEllipseDetectFromAll::~CEllipseDetectFromAll()
{

}

int CEllipseDetectFromAll::Detect(IplImage* pInput, vector<stEllipse>& ellipseArray, double fitErr, double minCount)
{
	//1. detect contours
	CvMemStorage* stor = cvCreateMemStorage(0);
	CvSeq* cont = cvCreateSeq(CV_SEQ_ELTYPE_POINT, sizeof(CvSeq), sizeof(CvPoint) , stor);
	int Nc = cvFindContours(pInput, stor, &cont);

	//2. remove the too small and long contours
	vector<POINT2> remainPts;
	CvPoint* PointArray;	

	IplImage* pDisp = cvCreateImage( cvGetSize(pInput), 8, 3 );
	memset(pDisp->imageData, 0, pDisp->height*pDisp->widthStep);

	for(; cont; cont = cont->h_next)
	{   
		printf(".");

		int i; 
		int count = cont->total; // This is number point in contour
		CvPoint center;
		CvSize size;

		// Number point must be more than or equal to 6 
		if( count < 16 )
			continue;

		PointArray = (CvPoint*)malloc( count*sizeof(CvPoint) );
		cvCvtSeqToArray(cont, PointArray, CV_WHOLE_SEQ);
		vector<POINT2> pts;
		for(i=0; i<count; i++)
		{
			POINT2 tp;
			tp.x = (float)PointArray[i].x;
			tp.y = (float)PointArray[i].y;
			pts.push_back(tp);
		}
		free(PointArray);

		double longestDis = CalcLongestDis(pts);
		if(fabs(count-longestDis)<(0.3*longestDis) )
			continue;

		int r = (double)( rand() ) / (double)(RAND_MAX) * 255;
		int g = (double)( rand() ) / (double)(RAND_MAX) * 255;
		int b = (double)( rand() ) / (double)(RAND_MAX) * 255;
		printf("%d %d %d \n", r,g,b);
		cvDrawContours( pDisp, cont, CV_RGB(r,g,b), CV_RGB(r,g,b),0,1,8);

		for(int i=0; i<pts.size(); i++)
		{
			remainPts.push_back(pts[i]);
		}
	}
	cvSaveImage("d:\\contours.jpg", pDisp);

	//3. random hough ellipse detection from all remained contours
	RandomHoughEllipse(remainPts, ellipseArray, fitErr, minCount);	

	return 1;
}