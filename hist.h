
#ifndef  COREDLL_HIST_H 
#define  COREDLL_HIST_H

//#include "cv.h"
//#include "highgui.h"


#define DLL_EXPORT  _declspec(dllexport)

#include "gdal_priv.h"

DLL_EXPORT  void   AutoLevelImage1(unsigned char* image, int ht, int wd);
//void DLL_ EXPORT HistMatchingColor(IplImage* src, IplImage* dst);
DLL_EXPORT  int GeoMosaicWithBalance(char** filesPath, int nFile, char* outpath, GDALDataType dataType = GDT_Byte);






#endif