#include "dodging.h"



/*  image dodging based on MASK, written by xiedonghai, 2013.8.9
*/
void MaskDodging(IplImage* pSrc, IplImage* pDst)
{
	int nChannel = pSrc->nChannels;

	if(nChannel==1)
	{
		/*
		//convert from byte to float type
		IplImage *pFloatImage = cvCreateImage( cvGetSize(pSrc), IPL_DEPTH_32F,1);
		cvConvertScale(pSrc, pFloatImage);		
		//gaussian filter
		IplImage* pSmooth = cvCloneImage(pFloatImage);
		cvSmooth(pFloatImage, pSmooth);	    
		//get the high frequency part of image
		IplImage* pDiff = cvCloneImage(pFloatImage);
		cvSub(pFloatImage, pSmooth, pDiff);		
		//cvSaveImage("d:\\highFrequecyImage.bmp", pDiff);
		*/

		IplImage* pSmooth = cvCloneImage(pSrc);
		cvSmooth(pSrc, pSmooth, 2, 81, 81,40);
		cvSaveImage("d:\\smooth.jpg", pSmooth);

		//IplImage* pDiff = cvCloneImage(pSrc);
		//cvSub(pSrc, pSmooth, pDiff);
		int ht = pDst->height;
		int wd = pDst->width;
		int scanwd = pDst->widthStep;
		for(int j=0; j<ht; j++)
			for(int i=0; i<wd; i++)
			{
				unsigned char originValue = (unsigned char)(pSrc->imageData[j*scanwd+i]);
				unsigned char smoothValue = (unsigned char)(pSmooth->imageData[j*scanwd+i]);
                pDst->imageData[j*scanwd+i] = 128 + (originValue-smoothValue);
			}

		//stretch

		/*
		//convert
		//cvConvertScale(pDiff, pDst);
		CvScalar bgMean;
		bgMean.val[0] = 128;
		bgMean.val[1] = 128;
		bgMean.val[2] = 128;
		bgMean.val[3] = 1;
		cvAddS(pDiff, bgMean, pDst);
		//cvSaveImage("d:\\dodge.bmp", pDst);
		*/
	}

	if(nChannel==3)
	{
		IplImage* pSmooth = cvCloneImage(pSrc);
		cvSmooth(pSrc, pSmooth, 2, 81, 81,40);
		cvSaveImage("d:\\smooth.jpg", pSmooth);

		int ht = pDst->height;
		int wd = pDst->width;
		int scanwd = pDst->widthStep;

		//calculate the mean value of image
		int meanR = 128;
		int meanG = 128;
		int meanB = 128;
		int i,j;

		int sumr = 0;
		int sumg = 0;
		int sumb = 0;
		for( j=0; j<ht; j++)
			for( i=0; i<wd; i++)
			{
				unsigned char or = (unsigned char)(pSrc->imageData[j*scanwd+i*3]);
				unsigned char og = (unsigned char)(pSrc->imageData[j*scanwd+i*3+1]);
				unsigned char ob = (unsigned char)(pSrc->imageData[j*scanwd+i*3+2]);			
				sumr += or;
				sumg += og;
				sumb += ob;
			}
		meanR = sumr/(ht*wd);
		meanG = sumg/(ht*wd);
		meanB = sumb/(ht*wd);

		//dodging
		for( j=0; j<ht; j++)
			for( i=0; i<wd; i++)
			{
				int or = (unsigned char)(pSrc->imageData[j*scanwd+i*3]);
				int sr = (unsigned char)(pSmooth->imageData[j*scanwd+i*3]);
				pDst->imageData[j*scanwd+i*3] = max( 0, min(255, meanR + (or-sr)) );

				int og = (unsigned char)(pSrc->imageData[j*scanwd+i*3+1]);
				int sg = (unsigned char)(pSmooth->imageData[j*scanwd+i*3+1]);
				pDst->imageData[j*scanwd+i*3+1] = max( 0, min(255, meanG + (og-sg)) );

				int ob = (unsigned char)(pSrc->imageData[j*scanwd+i*3+2]);
				int sb = (unsigned char)(pSmooth->imageData[j*scanwd+i*3+2]);
				pDst->imageData[j*scanwd+i*3+2] = max(0, min(255, meanB + (ob-sb)) );
			}
	}
}
