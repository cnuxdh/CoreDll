// CoreDll.cpp : Defines the entry point for the DLL application.
//

#include "stdafx.h"

#include "stdio.h"
#include "stdlib.h"

#include "../../Corelib/SVMCommon.h"

#include "exports.h"


#ifdef _MANAGED
#pragma managed(push, off)
#endif


//svm trained result
SVM_MODEL    model;
KERNEL_PARM  kernel_parm;

//weighted svm vector
int nFeat;
float* pSVMVec = NULL;


BOOL APIENTRY DllMain( HMODULE hModule,
                       DWORD  ul_reason_for_call,
                       LPVOID lpReserved
					 )
{
    return TRUE;
}

int LoadSVMResultWeight(char* svmTrainFile)
{
	FILE* fp = NULL;
	int nLen;
	long n;
	
	fp = fopen(svmTrainFile, "rb");
	fseek(fp,0,2);   /*将文件指针从文件头移动文件尾*/ 
	n=ftell(fp);       /*检测文件当前指针位置，求得文件长度*/ 
    fclose(fp);

    nLen = n / sizeof(float);

	if(pSVMVec != NULL) free(pSVMVec);

	nFeat = nLen;
	pSVMVec = (float*)malloc(sizeof(float)*nFeat);    
	fp = fopen(svmTrainFile, "rb");
	fread(pSVMVec, sizeof(float), nFeat, fp);
	fclose(fp);

	return 1;
}

double SVMClassifyWeight(float* pFeat, int ndim)
{
	double res;
	int i;
	res = 0;
	for(i=0; i<ndim; i++)
	{
		res += pFeat[i]*pSVMVec[i];
	}
	res += pSVMVec[ndim];

	return res;
}


int LoadSVMResult(char* svmTrainFile)
{
	ReadSVMTrainingRes(svmTrainFile, model, kernel_parm);
	return 1;
}

double SVMClassify(float* pFeat, int ndim)
{
	IN_SAM s;
	int i,j,k;
	double res;

	s.vector = (double*)malloc(sizeof(double)*ndim);

	for(j=0; j<ndim; j++)
		s.vector[j] = pFeat[j];		
	s.twonorm_sq = sprod_ss(s.vector, s.vector, ndim);
	res = classify_example(&model, &s, &kernel_parm, ndim);	

	free(s.vector);

	return res;
}


#ifdef _MANAGED
#pragma managed(pop)
#endif

