

#include "VideoWrite.h"


COpencvVideoWrite::COpencvVideoWrite()
{
	m_pWriter = NULL;
}

COpencvVideoWrite::~COpencvVideoWrite()
{

}

int COpencvVideoWrite::Create(char* filepath, int ht, int wd )
{
	int outCompressCodec = CV_FOURCC('X','V','I','D');
	int fps = 20;
	
	if(m_pWriter!=NULL)
		cvReleaseVideoWriter(&m_pWriter);
	
	m_pWriter=cvCreateVideoWriter(filepath, outCompressCodec, fps, cvSize(wd, ht), 1);
    
	if(m_pWriter == NULL)
		return 0;

	return 1;
}

int COpencvVideoWrite::WriteFrame(IplImage* pFrame)
{
	if(m_pWriter!=NULL)
		cvWriteFrame(m_pWriter, pFrame);
	
	return 1;
}

int COpencvVideoWrite::Close()
{
	if(m_pWriter!=NULL)
		cvReleaseVideoWriter(&m_pWriter);

	return 1;
}