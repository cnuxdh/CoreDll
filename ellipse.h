

#ifndef ELLIPSE_H
#define ELLIPSE_H

//corelib
#include "commondata.h"

//opencv
#include "cv.h"
#include "highgui.h"
#include "cvtypes.h"

#include "exports.h"


#include <vector>
using namespace std;


double CalcLongestDis(vector<POINT2> pts);


DLL_EXPORT void   DrawEllipse(stEllipse ellipsePara, IplImage* pImage);

DLL_EXPORT int    InvertEllipse(vector<POINT2> pts, stEllipse& ellipsePara, double& err);

DLL_EXPORT double EllipseFunc(double x, double y, stEllipse ellipsePara);

DLL_EXPORT int    DetectEllipse(IplImage* pInput, vector<stEllipse>& ellipseArray);

DLL_EXPORT int	  DetectEllipseUsingContour(IplImage* pInput, vector<stEllipse>& ellipseArray);

DLL_EXPORT double Pt2Ellipse(double x, double y, stEllipse ellipsePara);


class DLL_EXPORT CEllipseDetectBase
{
public:
	CEllipseDetectBase(){}
	virtual ~CEllipseDetectBase(){}

	virtual int Detect(IplImage* pInput, vector<stEllipse>& ellipseArray, double fitErr, double minCount){return 0;}
};

class DLL_EXPORT CEllipseDetectOneByOne: public CEllipseDetectBase
{
public:
	CEllipseDetectOneByOne();
	~CEllipseDetectOneByOne();
	
	int Detect(IplImage* pInput, vector<stEllipse>& ellipseArray, double fitErr, double minCount);
};

class DLL_EXPORT CEllipseDetectFromAll: public CEllipseDetectBase
{
public:
	CEllipseDetectFromAll();
	~CEllipseDetectFromAll();

	int Detect(IplImage* pInput, vector<stEllipse>& ellipseArray, double fitErr, double minCount);
};



#endif