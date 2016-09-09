
#include"windows.h"
#include "tchar.h"
#include <assert.h>
#include <time.h>
#include <string.h>

//#include "../../BundlerLib/BundlerLib/keys2a.h"
#include "keys2a.h"
#include "Matching.h"

#include "CommonFuncs1.h"

#include "../../CoreLib/ImageFunc.h"
#include "../../CoreLib/CommonFuncs.h"



void MathcingUsingSiftKey(int num_images, unsigned char** keys, int* num_keys, char* resFile)
{
	double ratio = 0.6;
	clock_t start;
	clock_t end ;    
	FILE* fp = NULL;
	int i;

	fp = fopen(resFile, "w");

	for( i = 0; i < num_images; i++) 
	{
		if (num_keys[i] == 0)
			continue;

		printf("[KeyMatchFull] Matching to image %d\n", i);

		start = clock();

		/* Create a tree from the keys */
		ANNkd_tree *tree = CreateSearchTree(num_keys[i], keys[i]);

		for (int j = 0; j < i; j++) 
		{
			if (num_keys[j] == 0)
				continue;

			/* Compute likely matches between two sets of keypoints */
			std::vector<KeypointMatch> matches = 
				MatchKeys(num_keys[j], keys[j], tree, ratio);

			int num_matches = (int) matches.size();

			if (num_matches >= 16) 
			{
				/* Write the pair */
				fprintf(fp, "%d %d\n", j, i);

				/* Write the number of matches */
				fprintf(fp, "%d\n", (int) matches.size());

				for (int i = 0; i < num_matches; i++) 
				{
					fprintf(fp, "%d %d\n", 
						matches[i].m_idx1, matches[i].m_idx2);
				}
			}
		}
		end = clock();    
		printf("[KeyMatchFull] Matching took %0.3fs\n", 
			(end - start) / ((double) CLOCKS_PER_SEC));
 	}
	fclose(fp);
}


int  ReadSiftKeys(FILE *fp, unsigned char **keys)
{
	int i, num, len;
	int j;

	std::vector<Keypoint *> kps;

	if (fscanf(fp, "%d %d", &num, &len) != 2)
	{
		printf("Invalid keypoint file\n");
		return 0;
	}

	if (len != 128) 
	{
		printf("Keypoint descriptor length invalid (should be 128).");
		return 0;
	}

	*keys = new unsigned char[128 * num + 8];

	unsigned char *p = *keys;
	for (i = 0; i < num; i++) 
	{
		float x, y, scale, ori;
		if (fscanf(fp, "%f %f %f %f\n", &y, &x, &scale, &ori) != 4) 
		{
			printf("Invalid keypoint file format.");
			return 0;
		}
		
		for(j=0; j<128; j++)
			fscanf(fp, "%hhu", p+j);
		
		p += 128;
	}

	return num; // kps;
}




/*
   binary image difference
   xoff,yoff: offset of the second image
*/
double TemplateMatching( IplImage* pL, IplImage* pR, int xoff, int yoff )
{
	int i,j;
	int ht = pL->height;
	int wd = pL->width;
	int scanwd = pL->widthStep;
	int lv,rv;

	double dif = 0;
	for(j=0; j<ht; j++)
		for(i=0; i<wd; i++)
		{
			//delete the pixel out of the circle
			double r = sqrt( (double)(j*j+i*i) );
			if( r>(wd*0.5) )
				continue;

			if( (j+yoff)<0 || (j+yoff)>=ht || (i+xoff)<0 || (i+xoff)>=wd)
				continue;

			int index1 = j*scanwd + i;
			int index2 = (j+yoff)*scanwd + i+xoff;
			lv = (unsigned char)(pL->imageData[index1]);
			rv = (unsigned char)(pR->imageData[index2]);
			
			dif += abs(lv-rv)/255;
		}

	return dif;
}


double LBPSimilarity(IplImage* src, IplImage* dst)
{
	unsigned char* pbuffer = NULL;
	int ht,wd;
	IplImageToGrayImage(src, &pbuffer, &ht, &wd);
	double lbpHist1[256];
	GenerateLBP(pbuffer, ht, wd, lbpHist1);
	free(pbuffer);

	//resize
	IplImageToGrayImage(dst, &pbuffer, &ht, &wd);
	double lbpHist2[256];
	GenerateLBP(pbuffer, ht, wd, lbpHist2);
	free(pbuffer);

	double sim = HistSimilar(lbpHist1, lbpHist2, 256, 0);

	return 1-sim;
}

/*
    similarity of two image with same size based on binary template matching
*/
double ImageSimilarity(IplImage* src, IplImage* dst)
{
	//segment
	cvAdaptiveThreshold(src, src, 255,
		CV_ADAPTIVE_THRESH_MEAN_C, 
		CV_THRESH_BINARY, 9, -5);
	cvAdaptiveThreshold(dst, dst, 255,
		CV_ADAPTIVE_THRESH_MEAN_C, 
		CV_THRESH_BINARY, 9, -5);       

	//similarity based on particle filter
	srand(1);
	double mindif = 10000000;
	for(int i=0; i<100; i++)
	{
		int xoff = ( (double)(rand())/(double)(RAND_MAX) - 0.5 )*8;
		int yoff = ( (double)(rand())/(double)(RAND_MAX) - 0.5 )*8;
		//printf("%d_%d  ", xoff, yoff);
		double dif = TemplateMatching(src, dst, xoff, yoff);  
		if(mindif>dif)
			mindif = dif;
	}
	return mindif;
}
