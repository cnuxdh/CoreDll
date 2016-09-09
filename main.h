
#ifndef COREDLL_MAIN_H
#define COREDLL_MAIN_H

#include "defs.h"

#include <vector>
using namespace std;


#define DLL_EXPORT  _declspec(dllexport)

int DLL_EXPORT siftFeatures(char* filename, struct feature** feat );
int DLL_EXPORT siftFeatures(char* filename, int dstHt, int dstWd, struct feature** feat );

void DLL_EXPORT PairMatchUsingSiftKey( unsigned char* lKey, int nLeftKey,
									  unsigned char* rKey, int nRightKey,
									  vector<int>& lKeyIds, vector<int>& rKeyIds);


#endif