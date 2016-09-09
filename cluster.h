#ifndef CLUSTER_H
#define CLUSTER_H

#include "exports.h"


#include <vector>
using namespace std;


//the pattern for cluster
typedef struct STRUPATTERN
{
	vector<double> s;
}stPattern;


//control parameters for isodata
typedef struct STRUISODATAPARA
{	
	int expectClusters;				//expected clustered number
	int initClusters;					//initial clustered number
	int minimalSamples;				//minimal number in one cluster
	double limitForSplit;			//the upper limit for split
	double limitForMerge;			//the low limit for merge
	double maximalNumberOfMerge;	//maximal number for merge
	int    maxIterations;			//the maximal iteration number
}stIsoDataParas;


//cluster 
typedef struct STRU_CLUSTER
{
	stPattern   centroid;             //the center patter
	double      std;
	stPattern   stdEach;              //std of each feature
	vector<int> clusterSet;         //sample index for the cluster            
}stCluster;


void DLL_EXPORT IsoData( vector<stPattern>& samples, vector<stCluster>& clusterRes, stIsoDataParas controlParas );


#endif