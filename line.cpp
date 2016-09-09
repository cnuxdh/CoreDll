

#include "line.h"


//corelib
#include "Corelib/CommonFuncs.h"


int HoughLines(IplImage* pImage, vector<stLINE>& lineArray)
{
	CvMemStorage* storage = cvCreateMemStorage(0);
	CvSeq* lines = 0;

	cvClearMemStorage(storage);
	lines = cvHoughLines2( pImage, storage, CV_HOUGH_PROBABILISTIC, 
		1,           //rho
		CV_PI/180,   //theta
		10,          //
		8,           //line length
		12);         //line gap
     
	for(int i=0; i<lines->total; i++)
	{
		CvPoint* line = (CvPoint*)cvGetSeqElem(lines,i);
		
		stLINE l;
		l.v1.x = line[0].x;
		l.v1.y = line[0].y;
		l.v2.x = line[1].x;
		l.v2.y = line[1].y;
		CalculateLinePolarParams(&l);
	
		lineArray.push_back(l);
	}
	
	cvReleaseMemStorage( &storage );

	return 1;
}

int HoughLines(IplImage* pImage, vector<stLINE>& lineArray,
			   double rho, double theta, int threshold,
			   double minLineLength, double maxLineGap)
{
	CvMemStorage* storage = cvCreateMemStorage(0);
	CvSeq* lines = 0;

	cvClearMemStorage(storage);
	lines = cvHoughLines2( pImage, storage, CV_HOUGH_PROBABILISTIC, 
		rho,         //rho
		theta,       //theta
		threshold,      //
		minLineLength,  //line 
		maxLineGap);    //line gap

	for(int i=0; i<lines->total; i++)
	{
		CvPoint* line = (CvPoint*)cvGetSeqElem(lines,i);

		stLINE l;
		l.v1.x = line[0].x;
		l.v1.y = line[0].y;
		l.v2.x = line[1].x;
		l.v2.y = line[1].y;
		CalculateLinePolarParams(&l);

		lineArray.push_back(l);
	}

	cvReleaseMemStorage( &storage );

	return 1;

}