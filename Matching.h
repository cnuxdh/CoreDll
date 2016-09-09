
//----------------------------------------------------------------------
// File:			Matching.h
// Programmer:		Xie Donghai
// Last modified:	06/28/2012 (Release 1.0)
// Description:		Basic include file for matching using approximate nearest
//					neighbor searching.
//----------------------------------------------------------------------
// Copyright (c)    All Rights Reserved.
//
//----------------------------------------------------------------------


#ifndef COREDLL_MATCHING_H
#define COREDLL_MATCHING_H

#include "exports.h"

#include "cv.h"
#include "highgui.h"
#include "cxcore.h"

/* Matching between images using sift keys based on ANN, 
   and matching result is saved into "resFile".
   All features of images are saved into "keys"      
input:
	num_images: the number of images
	keys:       all features of all images, each row represents features of one image
    num_keys:   the number of feature points of each image
	resFile:    file to save matching result
*/
void DLL_EXPORT MathcingUsingSiftKey(int num_images, unsigned char** keys, int* num_keys, char* resFile);

/* read sift key file format designed by Xie Donghai
input:
	fp:   key file
	keys: feature buffer
return: 
    the number of feature points 
*/
int  DLL_EXPORT ReadSiftKeys(FILE *fp, unsigned char **keys);



DLL_EXPORT double TemplateMatching( IplImage* pL, IplImage* pR, int xoff, int yoff );

DLL_EXPORT double ImageSimilarity(IplImage* src, IplImage* dst);

DLL_EXPORT double LBPSimilarity(IplImage* src, IplImage* dst);


#endif