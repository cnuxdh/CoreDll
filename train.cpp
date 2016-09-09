
#include "exports.h"


#include "../../Corelib/Gabor.h"
#include "../../Corelib/commonfile.h"
#include "../../Corelib/image.h"
#include "../../Corelib/Train.h"
#include "../../Corelib/ImageFunc.h"

#include "cv.h"
#include "highgui.h"
#include "cxcore.h"

#pragma comment (lib,"cv")
#pragma comment (lib,"highgui")
#pragma comment (lib,"cxcore")

/* Function: load image based on Opencv, and convert the image into 8-bit grey image   
*/	
void LoadImageUsingCV(char* filename, unsigned char** pImage, int* ht, int* wd)
{
	IplImage* pCvImage = NULL;
	int iht, iwd;
	int scanwd;
	int i,j;

	pCvImage = cvLoadImage(filename, 0);
	iht = pCvImage->height;
	iwd = pCvImage->width;
	scanwd = pCvImage->widthStep;
	
	*pImage = (unsigned char*)malloc(iht*iwd*sizeof(unsigned char));
	*ht = iht;
	*wd = iwd;

	for(j=0; j<iht; j++)
		for(i=0; i<iwd; i++)
		{
			(*pImage)[j*iwd+i] = pCvImage->imageData[j*scanwd+i];
		}

	cvReleaseImage(&pCvImage);
}


/* Function: generate feature samples from normalized images in input path
input:
	samPath:  directory including sample images and the size of each image is same
	featPath: directory to save feature file
	maskfile: mask file to select image point to generate feature
output:
	binary format *.dat files:  each ".dat" file save the feature vector generated from each image
return:
	return the dimension of feature vector
*/
int GenerateSampleFromImage(char* samPath,	char* featPath, char* maskfile)
{
	int nFeatDim = 0;
	int i,j,k;
	int jj,ii;
	char** filenames=NULL;
	int n,nfile;
	IplImage* pImage;	
	int** image;
	int ht, wd;
	int scanwd;
	float* pGaborFeat = NULL;
	float* pOutFeat = NULL;
	int  nGaborFeatDim;
	char featfile[256];
	FILE* fp = NULL;
	char* pdes;
	int index;
	unsigned char* pMask=NULL;
	GABOR_PARAM * pGaborParam = NULL;
	int nPtForTrain = 0;  //the point number which value is greater than 0 in the mask

	nGaborFeatDim = nscale*norient*recog_sam_ht*recog_sam_wd;

	pGaborParam = (GABOR_PARAM *)malloc(sizeof(GABOR_PARAM)*nscale);
	CalGaborParam(pGaborParam);
	
	//load mask file
	if(maskfile!=NULL)
	{
		LoadImageUsingCV(maskfile, &pMask, &ht, &wd);
		nPtForTrain = CalculateWhitePtNumber(pMask, ht, wd);
		nFeatDim = 2*nPtForTrain*nscale*norient;

		assert(ht==recog_sam_ht);
		assert(wd==recog_sam_wd);
	}
	else
	{
		nFeatDim = nGaborFeatDim;
	}

	//calculate the dim of Gabor feature containing all scales and orientations
	//nGaborFeatDim = nscale*norient*recog_sam_ht*recog_sam_wd;
	pGaborFeat = (float*)malloc(nGaborFeatDim*sizeof(float));
	pOutFeat = (float*)malloc(nFeatDim*sizeof(float));
   
	
	n = 0;
	nfile = 0;
	GetDirFileName(filenames, samPath, &n, &nfile, "jpg", 0);
	filenames = f2c(nfile, 512);
	GetDirFileName(filenames, samPath, &n, &nfile, "jpg", 1);

	for(i=0; i<nfile; i++)
	{
		printf("%s \n", filenames[i]);
		pImage = cvLoadImage(filenames[i], 0);

		if(pImage==NULL)
			continue;

		//image resize
		IplImage* pImageResize = cvCreateImage(cvSize(recog_sam_ht, recog_sam_ht), 8, 1);
		cvResize(pImage, pImageResize);	

		ht = pImageResize->height;
		wd = pImageResize->width;
		scanwd = pImageResize->widthStep;
		image = f2i(ht, wd);

		for(jj=0; jj<ht; jj++)
			for(ii=0; ii<wd; ii++)
			{
				image[jj][ii] = (unsigned char)(pImageResize->imageData[jj*scanwd+ii]);
			}
		//SaveBmp("d:\\sample.bmp", image, ht, wd);
		//image gray normalization
		//ProIlluminationFaceImage(image, ht, wd, NULL, NULL);
		//SaveBmp("d:\\sample_norm1.bmp", image, ht, wd);		

		//generate Gabor features
		FFTCalculateGaborFeat_Real(image, NULL, ht, wd, pGaborFeat, pGaborParam);
	
		if(pMask!=NULL)
		{
			CalAvgGaborRealFeat(pGaborFeat, pOutFeat, pMask);
			memcpy(pGaborFeat, pOutFeat, sizeof(float)*nFeatDim);
		}

		//save feature into binary file
		char title[512];
		char title1[512];
		strcpy(title, filenames[i]);
		pdes = strrchr(title, '\\');
		index = pdes - title + 1;
		strcpy(title1, title+index);		

		sprintf(featfile, "%s\\%s.dat", featPath, title1);
		fp = fopen(featfile, "wb");
		fwrite(pGaborFeat, sizeof(float), nFeatDim, fp);
		fclose(fp);

		FreeArray_int(image, ht, wd);
		cvReleaseImage(&pImage);
		cvReleaseImage(&pImageResize);
	}

	free(pGaborFeat);
	free(pOutFeat);
	
	FreeGaborParam(pGaborParam);

	return nFeatDim;
}


