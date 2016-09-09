
#include "stdio.h"
#include "image.h"

#include "gdal_priv.h"


//////////////////////////////////////////////////////////////////////////
CImageX::CImageX()
{
	m_pGrey = NULL;
	m_pCvImage = NULL;
	m_ht = 0;
	m_wd = 0;
}
CImageX::~CImageX()
{
	if(m_pCvImage!=NULL)
		cvReleaseImage(&m_pCvImage);
	if(m_pGrey!=NULL)
		free(m_pGrey);
}
int CImageX::Load(char* filepath)
{
	int i,j;

	if(m_pCvImage!=NULL)
		cvReleaseImage(&m_pCvImage);
	if(m_pGrey!=NULL)
		free(m_pGrey);
	
	strcpy(m_cFilePath, filepath);

	//load using opencv
	m_pCvImage = cvLoadImage(filepath, 0);
	if(m_pCvImage==NULL)
	{
		printf("Image can not be opened! \n");
		return 0;
	}

    m_ht = m_pCvImage->height;
	m_wd = m_pCvImage->width;

	//save as image matrix
	if(m_pGrey!=NULL)
		free(m_pGrey);
	m_pGrey = (unsigned char*)malloc(m_ht*m_wd);
	
	for(j=0; j<m_ht; j++)
		for(i=0; i<m_wd; i++)
		{
			m_pGrey[j*m_wd+i] = m_pCvImage->imageData[j*m_pCvImage->widthStep+i];
		}

	return 1;
}
int CImageX::GetHt()
{
	return m_ht;
}
int CImageX::GetWd()
{
	return m_wd;
}
//return the buffer matrix pointer
unsigned char* CImageX::GetBuffer()
{
	return m_pGrey;
}
//return the interface of opencv
IplImage* CImageX::GetCvArr()
{
	return m_pCvImage;
}
char* CImageX::GetFilePath()
{
	return m_cFilePath;
}
//////////////////////////////////////////////////////////////////////////





//////////////////////////////////////////////////////////////////////////
CImageXColor::CImageXColor()
{
	m_pGrey = NULL;
	m_pCvImageColor = NULL;
	m_pCvImageGray = NULL;
	m_ht = 0;
	m_wd = 0;
}
CImageXColor::~CImageXColor()
{
	if(m_pCvImageColor!=NULL)
		cvReleaseImage(&m_pCvImageColor);

	if(m_pCvImageGray!=NULL)
		cvReleaseImage(&m_pCvImageGray);

	if(m_pGrey!=NULL)
		free(m_pGrey);
}
int CImageXColor::Load(char* filepath)
{
	int i,j;

	if(m_pCvImageColor!=NULL)
		cvReleaseImage(&m_pCvImageColor);

	if(m_pCvImageGray!=NULL)
		cvReleaseImage(&m_pCvImageGray);

	if(m_pGrey!=NULL)
		free(m_pGrey);

	strcpy(m_cFilePath, filepath);

	//load using opencv
	m_pCvImageColor = cvLoadImage(filepath, 1);
	if(m_pCvImageColor==NULL)
	{
		printf("Image can not be opened! \n");
		return 0;
	}

	m_ht = m_pCvImageColor->height;
	m_wd = m_pCvImageColor->width;

	//save as image matrix
	if(m_pGrey!=NULL)
		free(m_pGrey);
	m_pGrey = (unsigned char*)malloc(m_ht*m_wd);

	m_pCvImageGray = cvLoadImage(filepath, 0);

	for(j=0; j<m_ht; j++)
		for(i=0; i<m_wd; i++)
		{
			m_pGrey[j*m_wd+i] = m_pCvImageGray->imageData[j*m_pCvImageGray->widthStep+i];
		}

	return 1;
}
int CImageXColor::GetHt()
{
	return m_ht;
}
int CImageXColor::GetWd()
{
	return m_wd;
}
//return the buffer matrix pointer
unsigned char* CImageXColor::GetBuffer()
{
	return m_pGrey;
}
//return the interface of opencv
IplImage* CImageXColor::GetCvArr()
{
	return m_pCvImageGray;
}
IplImage* CImageXColor::GetCvArrColor()
{
	return m_pCvImageColor;
}
char* CImageXColor::GetFilePath()
{
	return m_cFilePath;
}
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
CImageGdal::CImageGdal()
{
	mpCvImageGray = NULL;
	mpCvImageColor = NULL;
	mpGrey = NULL;
}
CImageGdal::~CImageGdal()
{
	//clear
	if(mpGrey!=NULL)
		free(mpGrey);
	if(mpCvImageGray!=NULL)
		cvReleaseImage(&mpCvImageGray);
	if(mpCvImageColor!=NULL)
		cvReleaseImage(&mpCvImageColor);

	for(int i=0; i<mBandData.size(); i++)
	{
		cvReleaseMat( &(mBandData[i]) );
	}
	mBandData.clear();
}

int CImageGdal::Load(char* filepath)
{
	//clear
	if(mpGrey!=NULL)
		free(mpGrey);
	if(mpCvImageGray!=NULL)
		cvReleaseImage(&mpCvImageGray);
	if(mpCvImageColor!=NULL)
		cvReleaseImage(&mpCvImageColor);

	for(int i=0; i<mBandData.size(); i++)
	{
		cvReleaseMat( &(mBandData[i]) );
	}
	mBandData.clear();

	char cSuffix[256];
	//retrieve the suffix
	char* pdes = strrchr(filepath, '.');
    int   index = pdes - filepath+1;
    strcpy(cSuffix, filepath+index);
    
	GDALAllRegister();

	GDALDriver* poDriver = NULL;
	GDALDataset *poDataset = NULL;   //GDAL数据集
	GDALRasterBand *poBand = NULL;

	char szFormat[256];

	if( strcmp(cSuffix, "tif")==0 )
	{
		strcpy(szFormat, "GTiff");
	}
	else if( strcmp(cSuffix, "jpg")==0 || strcmp(cSuffix, "JPG")==0 || strcmp(cSuffix, "jpeg")==0)
	{
		strcpy(szFormat, "JPEG");
	}
	else if( strcmp(cSuffix, "img")==0 )
	{
		strcpy(szFormat, "HFA");
	}	
	else
	{

	}

	poDriver = GetGDALDriverManager()->GetDriverByName(szFormat);
	if( poDriver == NULL )
	{
		printf("Gdal Open file failed ! \n");
		return 0;
	}

	poDataset = (GDALDataset *)GDALOpen(filepath, GA_ReadOnly );
	int nband = poDataset->GetRasterCount();
	printf("band number: %d \n", nband);

	poBand  = poDataset->GetRasterBand(1);
	int wd = poBand->GetXSize();
	int ht = poBand->GetYSize();
	GDALDataType nType = poBand->GetRasterDataType();
	const char* datatypename=GDALGetDataTypeName(nType);
	printf("Data type: %s \n", datatypename);
    
	switch (nType)
	{
	case GDT_Byte:
			for(int i=0; i<nband; i++)
			{
				poBand  = poDataset->GetRasterBand(i+1);
				int wd = poBand->GetXSize();
				int ht = poBand->GetYSize();
				unsigned char* pBuffer = (unsigned char*)malloc(ht*wd);		
				poBand->RasterIO(GF_Read,0,0,wd,ht,pBuffer,wd,ht,GDT_Byte,0,0);

				//save into mat
				//CvMat* pMat = cvCreateMat(ht, wd, CV_32FC1);
				CvMat* pMat = cvCreateMat(ht, wd, CV_8UC1);	
				if(pMat!=NULL)
				{
					//assert(pMat!=NULL);
					for(int m=0; m<ht; m++)
						for(int n=0; n<wd; n++)
						{
							//cvmSet(pMat, m, n, (double)(pBuffer[m*wd+n]) );
							cvSetReal2D(pMat, m, n, (double)(pBuffer[m*wd+n]));
						}
					mBandData.push_back(pMat);
				}
				free(pBuffer);
			}
			break;		
	case GDT_Int16:
			for(int i=0; i<nband; i++)
			{
				poBand  = poDataset->GetRasterBand(i+1);
				int wd = poBand->GetXSize();
				int ht = poBand->GetYSize();
				short* pBuffer = (short*)malloc(ht*wd*sizeof(short));		
				poBand->RasterIO(GF_Read,0,0,wd,ht,pBuffer,wd,ht,GDT_Int16,0,0);

				//save into mat
				CvMat* pMat = cvCreateMat(ht, wd, CV_32FC1);
				//CvMat* pMat = cvCreateMat(ht, wd, CV_8UC1);	
				if(pMat!=NULL)
				{
					//assert(pMat!=NULL);
					for(int m=0; m<ht; m++)
						for(int n=0; n<wd; n++)
						{
							//cvmSet(pMat, m, n, (double)(pBuffer[m*wd+n]) );
							cvSetReal2D(pMat, m, n, (double)(pBuffer[m*wd+n]));
						}
						mBandData.push_back(pMat);
				}

				free(pBuffer);
			}		
		break;
	case GDT_Float32:

		break;
	case GDT_Float64:

		break;
	default:
		break;
	}

	//transform into OpenCV format
	mpCvImageGray  = cvCreateImage(cvSize(wd, ht), 8, 1);
	mpCvImageColor = cvCreateImage(cvSize(wd, ht), 8, 3);	
	mpGrey = (unsigned char*)malloc(ht*wd);
	mHt = ht;
	mWd = wd;
	int scanwdGray = mpCvImageGray->widthStep;
	int scanwdColor = mpCvImageColor->widthStep;
	
	if(mBandData.size()==1)
	{		
		for(int j=0; j<ht; j++)
			for(int i=0; i<wd; i++)
			{
				mpCvImageGray->imageData[j*scanwdGray+i] = cvGetReal2D(mBandData[0], j, i);
				mpGrey[j*wd+i] = cvGetReal2D(mBandData[0], j, i);
				mpCvImageColor->imageData[j*scanwdColor+3*i]   = cvGetReal2D(mBandData[0], j, i);
				mpCvImageColor->imageData[j*scanwdColor+3*i+1] = cvGetReal2D(mBandData[0], j, i);
				mpCvImageColor->imageData[j*scanwdColor+3*i+2] = cvGetReal2D(mBandData[0], j, i);
			}
	}
	else 
	{
		for(int j=0; j<ht; j++)
			for(int i=0; i<wd; i++)
			{
				mpCvImageGray->imageData[j*scanwdGray+i] = cvGetReal2D(mBandData[0], j, i);
				mpGrey[j*wd+i] = cvGetReal2D(mBandData[0], j, i);
				mpCvImageColor->imageData[j*scanwdColor+3*i]   = cvGetReal2D(mBandData[2], j, i);
				mpCvImageColor->imageData[j*scanwdColor+3*i+1] = cvGetReal2D(mBandData[1], j, i);
				mpCvImageColor->imageData[j*scanwdColor+3*i+2] = cvGetReal2D(mBandData[0], j, i);
			}
	}
	
	GDALClose( (GDALDatasetH) poDataset );

	return 1;
}

int CImageGdal::GetHt()
{
	return mHt;
}
int CImageGdal::GetWd()
{
	return mWd;
}

char* CImageGdal::GetFilePath()
{
	return NULL;
}
unsigned char* CImageGdal::GetBuffer()
{
	return mpGrey;
}
IplImage* CImageGdal::GetCvArr()
{
	return mpCvImageGray;
}
IplImage* CImageGdal::GetCvArrColor()
{
	return mpCvImageColor;
}

int CImageGdal::GetBandCount()
{
	return mBandData.size();
}
//get the value of special channel
double CImageGdal::GetValue(int row, int col, int band)
{
	return cvGetReal2D( mBandData[band], row, col );
	//return 0;
}
//////////////////////////////////////////////////////////////////////////


//
int LoadGrayImageGeneral(char* filepath, double** pDst, int* dstHt, int* dstWd)
{
	int  i,j;
	char cSuffix[256];

	//retrieve the suffix
	char* pdes = strrchr(filepath, '.');
	int   index = pdes - filepath+1;
	strcpy(cSuffix, filepath+index);

	GDALAllRegister();

	GDALDriver* poDriver   = NULL;
	GDALDataset *poDataset = NULL;   //GDAL数据集
	GDALRasterBand *poBand = NULL;

	char szFormat[256];

	if( strcmp(cSuffix, "tif")==0 )
	{
		strcpy(szFormat, "GTiff");
	}
	else if( strcmp(cSuffix, "jpg")==0 || strcmp(cSuffix, "JPG")==0 || strcmp(cSuffix, "jpeg")==0)
	{
		strcpy(szFormat, "JPEG");
	}
	else if( strcmp(cSuffix, "img")==0 )
	{
		strcpy(szFormat, "HFA");
	}	
	else
	{

	}

	poDriver = GetGDALDriverManager()->GetDriverByName(szFormat);
	if( poDriver == NULL )
	{
		printf("Gdal Open file failed ! \n");
		return 0;
	}

	poDataset = (GDALDataset *)GDALOpen(filepath, GA_ReadOnly );
	int nband = poDataset->GetRasterCount();
	printf("band number: %d \n", nband);

	poBand  = poDataset->GetRasterBand(1);
	int wd = poBand->GetXSize();
	int ht = poBand->GetYSize();
	GDALDataType nType = poBand->GetRasterDataType();
	const char* datatypename=GDALGetDataTypeName(nType);
	printf("Data type: %s \n", datatypename);


	*dstHt = ht;
	*dstWd = wd;
	(*pDst) = (double*)malloc(ht*wd*sizeof(double));

	if((*pDst) == NULL)
	{
		printf("Memory is not enough! \n");
		return 0;
	}

	switch (nType)
	{
	case GDT_Byte:
		//for( i=0; i<nband; i++)
		i=0;
		{
			poBand  = poDataset->GetRasterBand(i+1);
			int wd = poBand->GetXSize();
			int ht = poBand->GetYSize();
			unsigned char* pBuffer = (unsigned char*)malloc(ht*wd);		
			poBand->RasterIO(GF_Read,0,0,wd,ht,pBuffer,wd,ht,GDT_Byte,0,0);
            
			for(i=0; i<wd*ht; i++)
				(*pDst)[i] = pBuffer[i];
			
			free(pBuffer);
		}
		break;		
	case GDT_Int16:
		//for( i=0; i<nband; i++)
		i = 0;
		{
			poBand  = poDataset->GetRasterBand(i+1);
			int wd = poBand->GetXSize();
			int ht = poBand->GetYSize();
			short* pBuffer = (short*)malloc(ht*wd*sizeof(short));		
			poBand->RasterIO(GF_Read,0,0,wd,ht,pBuffer,wd,ht,GDT_Int16,0,0);

			for(i=0; i<wd*ht; i++)
				(*pDst)[i] = pBuffer[i];

			free(pBuffer);
		}		
		break;
	case GDT_Float32:
		i = 0;
		{
			poBand  = poDataset->GetRasterBand(i+1);
			int wd = poBand->GetXSize();
			int ht = poBand->GetYSize();
			float* pBuffer = (float*)malloc(ht*wd*sizeof(float));		
			poBand->RasterIO(GF_Read, 0, 0, wd, ht, pBuffer, wd, ht, GDT_Float32, 0, 0);

			for(i=0; i<wd*ht; i++)
				(*pDst)[i] = pBuffer[i];

			free(pBuffer);
		}

		break;
	case GDT_Float64:

		break;
	default:
		break;
	}

	GDALClose( (GDALDatasetH) poDataset );

	return 1;
}


//using opencv
int SaveJpeg(char* filename, unsigned char* pbuffer, int ht, int wd)
{
	IplImage* pImage = cvCreateImage(cvSize(wd, ht), 8, 1);
    int scanwd = pImage->widthStep;
    
	for(int j=0; j<ht; j++)
		for(int i=0; i<wd; i++)
		{
			pImage->imageData[j*scanwd+i] = pbuffer[j*wd+i];
		}

	cvSaveImage(filename, pImage);
	cvReleaseImage(&pImage);

	return 1;
}