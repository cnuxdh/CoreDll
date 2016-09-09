
#include "corelib/commondata.h"
#include "corelib/commonfuncs.h"
#include "Corelib/image.h"

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "SuperResolution.h"

void POCS(double* x, double* y, int* v, int npt, char* outfile)
{
	double scale = 2.0;
	int    ht,wd;
	int    i,j;

	double minx = 1000000;
	double maxx = -1000000;
	double miny = 1000000;
	double maxy = -1000000;

	//determine the image size
	for(i=0; i<npt; i++)
	{
		if( minx>x[i] )	minx = x[i];
		if( maxx<x[i] )	maxx = x[i];
		if( miny>y[i] )	miny = y[i];
		if( maxy<y[i] )	maxy = y[i];
	}
	minx = max(0, minx);
	miny = max(0, miny);

	//get the image size of superresolution
	ht = int(maxy-miny+0.5)*scale;
	wd = int(maxx-minx+0.5)*scale;

	double* sv = (double*)malloc(ht*wd*sizeof(double));
	double* prev = (double*)malloc(ht*wd*sizeof(double));
	memset(sv, 0, ht*wd*sizeof(double));
	memset(prev, 0, ht*wd*sizeof(double));

	//filter kernel
	/*
	double kernel[25] = {   0.0119, 0,      0.0476, 0,      0.0119,
							0,      0.0476, 0.0952, 0.0476, 0,
							0.0476, 0.0952, 0.1905, 0.0952, 0.0476,
							0,      0.0476, 0.0952, 0.0476, 0,
							0.0119, 0,      0.0476, 0,      0.0119};
	int    kht = 5;
	int    kwd = 5;   
	*/
	
	//gaussian kernel, sigma=0.33	
	double kernel[9] = {    0.0001,    0.0097,   0.0001,
							0.0097,    0.9606,   0.0097,
							0.0001,    0.0097,   0.0001};
    int kht = 3;
	int kwd = 3;
							
	/*
	double kernel[25] = {   0.0001,    0.0020,    0.0055 ,   0.0020,    0.0001,
							0.0020 ,   0.0422,    0.1171 ,   0.0422,    0.0020,
							0.0055 ,  0.1171 ,   0.3248 ,   0.1171 ,   0.0055,
							0.0020 ,   0.0422 ,   0.1171,    0.0422 ,   0.0020,
							0.0001 ,   0.0020  ,  0.0055,    0.0020 ,   0.0001
						};
	int    kht = 5;
	int    kwd = 5;
	*/
	
	//initialize 
	for(i=0; i<npt; i++)
	{
		int tx = x[i]*scale+0.5;
		int ty = y[i]*scale+0.5;
		if(tx>=0 && tx<wd && ty>=0 && ty<ht)
		{
			sv[ty*wd+tx] = v[i];
		}
	}	
	SaveBmp("d:\\image.bmp", sv, ht, wd);
	memcpy( prev, sv, sizeof(double)*ht*wd ); 

	double error=0;
	int    iteration = 0;
	do 
	{
		iteration++;		

		Filter2D(sv, ht, wd, kernel, kht, kwd);		
		
		//SaveBmp("d:\\sr.bmp", sv, ht, wd);
		//SaveBmp("d:\\prev.bmp", prev, ht, wd);

		//error
		error=0;
		for(i=0; i<(ht*wd); i++)
		{
			error += fabs(sv[i]-prev[i]);
		}
		error /= (ht*wd);
		
		//
		for(i=0; i<npt; i++)
		{
			int tx = x[i]*scale+0.5;
			int ty = y[i]*scale+0.5;
			if(tx>=0 && tx<wd && ty>=0 && ty<ht)
			{
				sv[ty*wd+tx] = v[i];
			}
		}
		//save 
		memcpy( prev, sv, sizeof(double)*ht*wd );

		printf("POCS iteration: %d  error: %lf \n", iteration, error);

	}while( error>0.001 && iteration<50);

	printf("POCS finished! \n");
	Filter2D(sv, ht, wd, kernel, kht, kwd);	
	SaveBmp(outfile, sv, ht, wd);	

	free(sv);
	free(prev);
}
