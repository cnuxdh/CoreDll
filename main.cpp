
#include "main.h"

#include "vim_sift.h"
#include "keys2a.h"


int siftFeatures(char* filename, int dstHt, int dstWd, struct feature** feat )
{
	IplImage* img = cvLoadImage(filename);

	//resize the image
	IplImage* pSmallImage = cvCreateImage(cvSize(dstWd, dstHt), img->depth, img->nChannels);
	cvResize(img, pSmallImage);

	
	int res = _sift_features( pSmallImage, feat, SIFT_INTVLS, SIFT_SIGMA, SIFT_CONTR_THR,
		SIFT_CURV_THR, SIFT_IMG_DBL, SIFT_DESCR_WIDTH,
		SIFT_DESCR_HIST_BINS );		
	
	cvReleaseImage(&pSmallImage);
	cvReleaseImage(&img);

	return res;
}


int siftFeatures(char* filename, struct feature** feat )
{
	IplImage* img = cvLoadImage(filename);

	return _sift_features( img, feat, SIFT_INTVLS, SIFT_SIGMA, SIFT_CONTR_THR,
		SIFT_CURV_THR, SIFT_IMG_DBL, SIFT_DESCR_WIDTH,
		SIFT_DESCR_HIST_BINS );

	cvReleaseImage(&img);
}

void PairMatchUsingSiftKey(unsigned char* lKey, int nLeftKey,
						  unsigned char* rKey, int nRightKey,
						  vector<int>& lKeyIds, vector<int>& rKeyIds)
{
	double ratio = 0.7;

	/* Create a tree from the keys */
	//ANNkd_tree *tree = CreateSearchTree(nRightKey, rKey);	
	
	ANNpointArray pts = annAllocPts(nRightKey, 128);
	for (int i = 0; i < nRightKey; i++) 
	{
		memcpy(pts[i], rKey + 128 * i, sizeof(unsigned char) * 128);
	}
	
	/* Create a search tree for k2 */
	ANNkd_tree *tree = new ANNkd_tree(pts, nRightKey, 128, 16);

	/* Compute likely matches between two sets of keypoints */
	std::vector<KeypointMatch> matches = MatchKeys( nLeftKey, lKey, tree, ratio);

	int num_matches = (int)matches.size();
	for (int i = 0; i < num_matches; i++) 
	{
		//fprintf(fp, "%d %d\n", matches[i].m_idx1, matches[i].m_idx2);
		lKeyIds.push_back( matches[i].m_idx1);
		rKeyIds.push_back( matches[i].m_idx2);
	}

	annDeallocPts(pts);
	delete tree;
}
