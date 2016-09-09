

#ifndef EXPORTS_H
#define EXPORTS_H

	

#ifdef _WIN32
# define DLL_EXPORT __declspec( dllexport )
#else
# define DLL_EXPORT
#endif



DLL_EXPORT int    GenerateSampleFromImage(char* samPath, char* featPath, char* maskfile);
DLL_EXPORT void   LoadImageUsingCV(char* filename, unsigned char** pImage, int* ht, int* wd);
DLL_EXPORT int    LoadSVMResult(char* svmTrainFile);
DLL_EXPORT double SVMClassify(float* pFeat, int ndim);
DLL_EXPORT int    LoadSVMResultWeight(char* svmTrainFile);
DLL_EXPORT double SVMClassifyWeight(float* pFeat, int ndim);




#endif