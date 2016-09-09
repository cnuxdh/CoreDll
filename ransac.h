#ifndef RANSAC_H
#define RANSAC_H


#include "cv.h"
#include "highgui.h"


void RANSAC_homography(int num, CvPoint2D64f *m1, CvPoint2D64f *m2, CvMat *H,CvMat *inlier_mask);



#endif