
#ifndef TRAIN_H
#define TRAIN_H
	

//void GenerateMask(int ht, int wd, int nStep, char* maskfile);
//int  CalculateWhitePtNumber(unsigned char* pMask, int ht, int wd);
//void LoadImageUsingCV(char* filename, unsigned char** pImage, int* ht, int* wd);	
int  GenerateSampleFromImage(char* samPath, char* featPath, char* maskfile, GABOR_PARAM * pGaborParam);



#endif