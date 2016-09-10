#ifndef RANSAC_H
#define RANSAC_H

#ifdef OPENCV_1X 
	#include "cv.h"
	#include "highgui.h"
	#include "cxcore.h"
#else
	#include "opencv2/core/core.hpp"
	#include "opencv2/highgui/highgui.hpp"
	#include "opencv2/calib3d/calib3d.hpp"
	using namespace cv;
#endif
//#include "cv.h"
//#include "highgui.h"


void RANSAC_homography(int num, CvPoint2D64f *m1, CvPoint2D64f *m2, CvMat *H,CvMat *inlier_mask);



#endif