

#ifndef LINE_H
#define LINE_H


#include "exports.h"

//corelib
#include "commondata.h"

//opencv
#include "cv.h"
#include "highgui.h"
#include "cvtypes.h"

#include <vector>
using namespace std;

DLL_EXPORT int HoughLines(IplImage* pImage, vector<stLINE>& lineArray);

DLL_EXPORT int HoughLines(IplImage* pImage, vector<stLINE>& lineArray,
						  double rho, double theta, int threshold,
						  double minLineLength=0, double maxLineGap=0);


#endif