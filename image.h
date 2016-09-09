#ifndef MYIMAGE_H
#define MYIMAGE_H

#include "exports.h"

#include<vector>

#include "cv.h"
#include "highgui.h"
#include "cvtypes.h"

using namespace std;


class DLL_EXPORT CImageParent
//class  CImageParent
{
public:
	CImageParent(){}
	virtual ~CImageParent(){}

	virtual int Load(char* filepath){return 1;}
	virtual int GetHt(){return 0;}
	virtual int GetWd(){return 0;}
	virtual int GetBandCount(){return 0;}
	virtual char* GetFilePath(){return NULL;}

	//get the value of special channel
	virtual double GetValue(int row, int col, int band){return 0;}

	//get data buffer
	virtual unsigned char* GetBuffer(){return NULL;}
	virtual IplImage* GetCvArr(){return NULL;}      //return 8 bit gray image
	virtual IplImage* GetCvArrColor(){return NULL;} //return original image
};

//for 8 bits gray image processing
class DLL_EXPORT CImageX: public CImageParent
{
public:
	 CImageX();
	~CImageX();

	int Load(char* filepath);
	int GetHt();
	int GetWd();
	char* GetFilePath();
	unsigned char* GetBuffer();
	IplImage* GetCvArr();
	IplImage* GetCvArrColor(){return NULL;}

private:
	unsigned char*  m_pGrey;
	int        m_ht,m_wd;
	IplImage*  m_pCvImage;
	char m_cFilePath[256];
};

//for 8 bits and 24 bits image processing
class DLL_EXPORT CImageXColor: public CImageParent
{
public:
	CImageXColor();
	~CImageXColor();

	int Load(char* filepath);
	int GetHt();
	int GetWd();
	char* GetFilePath();
	unsigned char* GetBuffer();
	IplImage* GetCvArr();
	IplImage* GetCvArrColor();

private:
	unsigned char*  m_pGrey;
	int        m_ht,m_wd;
	IplImage*  m_pCvImageColor;
	IplImage*  m_pCvImageGray;
	char m_cFilePath[256];
};

//for multiple channels, different type, based on Gdal
class DLL_EXPORT CImageGdal: public CImageParent
{
public:
	CImageGdal();
	~CImageGdal();

	int Load(char* filepath);
	int GetHt();
	int GetWd();
	char* GetFilePath();

	unsigned char* GetBuffer();
	IplImage* GetCvArr();
	IplImage* GetCvArrColor();

	int GetBandCount();
	
	//get the value of special channel
	double GetValue(int row, int col, int band);

private:
	int mHt,mWd;  //the dimension of image	
	int mBand;    //the band number of image
	
	//for data 
	vector<CvMat*> mBandData;

	//for opencv
	IplImage*  mpCvImageColor;
	IplImage*  mpCvImageGray;

	//for bitmap
	unsigned char*  mpGrey;
};

//single function to read image generally, especially for huge image data, added by Donghai Xie, 2013.7.10
DLL_EXPORT int LoadGrayImageGeneral(char* filepath, double** pDst, int* dstHt, int* dstWd);
DLL_EXPORT int SaveJpeg(char* filename, unsigned char* pbuffer, int ht, int wd);



#endif