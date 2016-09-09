
#ifndef  COREDLL_HIST_H 
#define  COREDLL_HIST_H

//#include "cv.h"
//#include "highgui.h"


#define DLL_EXPORT  _declspec(dllexport)

#include "gdal_priv.h"


//void DLL_EXPORT HistMatchingColor(IplImage* src, IplImage* dst);
DLL_EXPORT  int GeoMosaicWithBalance(char** filesPath, int nFile, char* outpath, GDALDataType dataType = GDT_Byte);






#endif