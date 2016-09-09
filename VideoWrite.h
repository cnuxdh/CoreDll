
#ifndef VIDEO_WRITE_H
#define VIDEO_WRITE_H

#include "exports.h"

//#include "cv.h"
#include "highgui.h"

//interface
class DLL_EXPORT CVideoWrite
{
public:
	CVideoWrite(){}
	virtual ~CVideoWrite(){}

	virtual int  Create(char* filepath, int ht, int wd ){return 1;}
	virtual int  WriteFrame(IplImage* pFrame){return 1;}
	virtual int  Close(){return 1;}
private:	
};


//video write using opencv functions
class DLL_EXPORT COpencvVideoWrite: public CVideoWrite
{
public:
	COpencvVideoWrite();
	~COpencvVideoWrite();
	int  Create(char* filepath, int ht, int wd );
	int  WriteFrame(IplImage* pFrame);
	int  Close();
private:
	CvVideoWriter *m_pWriter;
};


#endif