

#ifndef EXPORTS_H
#define EXPORTS_H

	
#define DLL_EXPORT  _declspec(dllexport)


DLL_EXPORT int    GenerateSampleFromImage(char* samPath, char* featPath, char* maskfile);
DLL_EXPORT void   LoadImageUsingCV(char* filename, unsigned char** pImage, int* ht, int* wd);
DLL_EXPORT int    LoadSVMResult(char* svmTrainFile);
DLL_EXPORT double SVMClassify(float* pFeat, int ndim);
DLL_EXPORT int    LoadSVMResultWeight(char* svmTrainFile);
DLL_EXPORT double SVMClassifyWeight(float* pFeat, int ndim);




#endif