
#include "math.h"
#include "time.h"


//corelib
#include "Corelib/CommonFuncs.h"


#include "cluster.h"


//cluster sampleIndex among samples into clusterRes
void Cluster(vector<stPattern>& samples, vector<int>& sampleIndex, vector<stCluster>& clusterRes)
{
	for(int i=0; i<sampleIndex.size(); i++)
	{
		int si = sampleIndex[i];

		double minDis = 100000;
		int    index = 0;

		for(int j=0; j<clusterRes.size(); j++)
		{
			double dis = 0; //VecDistance( samples[si].s, clusterRes[j].centroid.s );
			if(dis<minDis)
			{
				minDis = dis;
				index = j;
			}				
		}		
		clusterRes[index].clusterSet.push_back(si);
	}
}

//calculate the centor, std of the cluster
//samples: all original samples
void CalculateClusterParas(vector<stPattern>& samples, vector<stCluster>& clusterRes)
{
	for(int i=0; i<clusterRes.size(); i++)
	{
		int nClusterSize = clusterRes[i].clusterSet.size(); 
		int dim = clusterRes[i].centroid.s.size();

		stPattern center;
		center.s.resize(dim, 0);
		for(int j=0; j<nClusterSize; j++)
		{
			int index = clusterRes[i].clusterSet[j];
			for(int k=0; k<dim; k++)
				center.s[k] += samples[index].s[k];
		}
		//center
		for(int k=0; k<dim; k++)
			center.s[k] /= (double)(clusterRes[i].clusterSet.size());
		
		//std: whole vector
		double std = 0;
		for(int j=0; j<nClusterSize; j++)
		{
			int index = clusterRes[i].clusterSet[j];
			double std1 = 0;
			for(int k=0; k<dim; k++)
				std1 += ( center.s[k] - samples[index].s[k] )*( center.s[k] - samples[index].s[k] );
			std += std1;
		}
		std = sqrt( std/(double)( nClusterSize ) );

		//std: each feature
        stPattern stdEach;
		stdEach.s.resize(dim, 0);
		for(int j=0; j<nClusterSize; j++)
		{
			int index = clusterRes[i].clusterSet[j];
			double std1 = 0;
			for(int k=0; k<dim; k++)
			{
				stdEach.s[k] += ( center.s[k] - samples[index].s[k] )*( center.s[k] - samples[index].s[k] );
			}			
		}
		for(int k=0; k<dim; k++)
		{
			stdEach.s[k] = sqrt( stdEach.s[k]/(double)( nClusterSize ) );;
		}

		clusterRes[i].centroid = center;
		clusterRes[i].std = std;
		clusterRes[i].stdEach = stdEach;
	}
}

void Split(vector<stPattern>& samples, vector<stCluster>& clusterRes, stIsoDataParas controlParas)
{
	vector<int> sampleIndex;

	//split
	for(int i=0; i<clusterRes.size(); i++)
	{
		if( clusterRes[i].std > controlParas.limitForSplit )
		{
			vector<int> curSamIndex;
			curSamIndex = clusterRes[i].clusterSet;

			stCluster c1,c2;
			int dim = clusterRes[i].centroid.s.size();
			c1.centroid.s.resize(dim, 0);
			c2.centroid.s.resize(dim, 0);
			for(int k=0; k<dim; k++)
			{
				c1.centroid.s[k] = clusterRes[i].centroid.s[k] + clusterRes[i].stdEach.s[k]*0.5;
				c2.centroid.s[k] = clusterRes[i].centroid.s[k] - clusterRes[i].stdEach.s[k]*0.5;
			}

			//erase the original cluster
			clusterRes.erase( clusterRes.begin()+i );

			//add new clusters
			clusterRes.push_back(c1);
			clusterRes.push_back(c2);

			Cluster(samples, curSamIndex, clusterRes);

			//
			CalculateClusterParas(samples, clusterRes);


			i--;
		}
	}

	for(int ki=0; ki<clusterRes.size(); ki++)
	{
		if( clusterRes[ki].clusterSet.size()<controlParas.minimalSamples  )
		{
			//remove current cluster
			sampleIndex = clusterRes[ki].clusterSet;
			clusterRes.erase( clusterRes.begin()+ki );

			//cluster again
			Cluster(samples, sampleIndex, clusterRes);
			ki--;
		}
	}
	CalculateClusterParas(samples, clusterRes);
}

void Merge(vector<stPattern>& samples, vector<stCluster>& clusterRes, stIsoDataParas controlParas)
{
	for(int i=0; i<clusterRes.size(); i++)
	{
		for(int j=i+1; j<clusterRes.size(); j++)
		{
			int dim = clusterRes[i].centroid.s.size();
			double dis = 0; //VecDistance( clusterRes[i].centroid.s,  clusterRes[j].centroid.s);
			if(dis < controlParas.limitForMerge)
			{
				for(int k=0; k<dim; k++)
				{
					clusterRes[i].centroid.s[k] = (clusterRes[i].centroid.s[k] + clusterRes[j].centroid.s[k]) * 0.5;
				}
				for(int k=0; k<clusterRes[j].clusterSet.size(); k++)
				{
					clusterRes[i].clusterSet.push_back( clusterRes[j].clusterSet[k] );
				}
			
				clusterRes.erase( clusterRes.begin() + j );
				j--;
			}			
		}
	}
	CalculateClusterParas(samples, clusterRes);
}


//isodata algorithm, written by Xie Donghai, 2014.5.26
void IsoData( vector<stPattern>& samples, vector<stCluster>& clusterRes, stIsoDataParas controlParas )
{

	srand( time(NULL) );
	
	int nSamples = samples.size();

	// genrate cluster centers random
    vector<int> randSamIndex;
    for(int i=0; i<controlParas.initClusters; i++)
	{
		int ri = (double)(rand()) / (double)( RAND_MAX ) * (nSamples-1);
		randSamIndex.push_back(ri);
	}

	//initial cluster
	for(int i=0; i<randSamIndex.size(); i++)
	{
		stCluster c;
		c.centroid = samples[ randSamIndex[i] ];
		clusterRes.push_back(c);
	}

	//initialzing clustering
	vector<int> sampleIndex;
	for(int i=0; i<samples.size(); i++)
		sampleIndex.push_back(i);
	Cluster(samples, sampleIndex, clusterRes);

	//make sure the size of each cluster is more than the threthod
	for(int i=0; i<clusterRes.size(); i++)
	{
		if( clusterRes[i].clusterSet.size()<controlParas.minimalSamples  )
		{
			//remove current cluster
			sampleIndex = clusterRes[i].clusterSet;
			clusterRes.erase( clusterRes.begin()+i );
			//cluster again
			Cluster(samples, sampleIndex, clusterRes);
			i--;
		}
	}

	//calculate the centroid, std
	CalculateClusterParas(samples, clusterRes);
   

	int iterationNum = 0;

	while( iterationNum<controlParas.maxIterations )
	{
		int nc = clusterRes.size();

		printf("Iteration: %d \n", iterationNum);
		for(int k=0; k<clusterRes.size(); k++)
		{
			for(int m=0; m<clusterRes[k].centroid.s.size(); m++)
				printf("%lf ", clusterRes[k].centroid.s[m]);
			printf("\n");
		}

		if( nc<= (controlParas.expectClusters*0.5) )
		{
			Split(samples, clusterRes, controlParas);
		}
		else if( nc < (controlParas.expectClusters*2) )
		{
			if( iterationNum%2==0 )
			{
				Split(samples, clusterRes, controlParas);
			}
			else
			{
				Merge(samples, clusterRes, controlParas);
			}
		}
		//merge
		else //( nc >= (controlParas.expectClusters*2) )
		{
			Merge(samples, clusterRes, controlParas);
		}
		iterationNum++;
	}

}




