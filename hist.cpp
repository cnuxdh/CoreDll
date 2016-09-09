#include "hist.h"

//#include "Corelib/ImageFunc.h"
#include "geotiff.h"

#include "gdal_priv.h"
#include "ogr_spatialref.h"


#include <string>
#include <vector>

using namespace std;


/* histogram matching for color image   
*/
/*
void HistMatchingColor(IplImage* src, IplImage* dst)
{
	IplImage* rr = cvCreateImage(cvGetSize(src), src->depth, 1);
	IplImage* rg = cvCreateImage(cvGetSize(src), src->depth, 1);
	IplImage* rb = cvCreateImage(cvGetSize(src), src->depth, 1);

	IplImage* tr = cvCreateImage(cvGetSize(dst), dst->depth, 1);
	IplImage* tg = cvCreateImage(cvGetSize(dst), dst->depth, 1);
	IplImage* tb = cvCreateImage(cvGetSize(dst), dst->depth, 1);

	cvSplit( src, rr, rg, rb, NULL);
	cvSplit( dst, tr, tg, tb, NULL );

	HistMatching( (unsigned char*)(rr->imageData), 
		rr->height, rr->width, rr->widthStep, 
		(unsigned char*)(tr->imageData), 
		tr->height, tr->width, tr->widthStep);

	HistMatching( (unsigned char*)(rg->imageData), 
		rg->height, rg->width, rg->widthStep, 
		(unsigned char*)(tg->imageData), 
		tg->height, tg->width, tg->widthStep);


	HistMatching( (unsigned char*)(rb->imageData), 
		rb->height, rb->width, rb->widthStep, 
		(unsigned char*)(tb->imageData), 
		tb->height, tb->width, tb->widthStep);

	cvMerge(tr,tg,tb,NULL,dst);

	cvReleaseImage(&rr);
	cvReleaseImage(&rg);
	cvReleaseImage(&rb);
	cvReleaseImage(&tr);
	cvReleaseImage(&tg);
	cvReleaseImage(&tb);
}
*/



/*功能：对图像进行auto level运算
默认黑白各自去掉 0.5% 的范围
*/
void   AutoLevelImage1(unsigned char* image, int ht, int wd)
{
	int i,j;
	int index;
	int hist[256];
	int imgSize;
	float delRatio = 0.01;
	int   delSum;
	int   s;
	int   lg,rg;
	float scale;

	memset(hist, 0, sizeof(int)*256);

	imgSize = ht*wd;
	int sum = 0;
	//1: 统计灰度分布的直方图
	for(i=0; i<imgSize; i++)
	{
		if(image[i]>0 )
		{
			hist[ image[i] ] ++;
			sum++;
		}
	}

	//2: 搜索最黑与最白的边界
	delSum = sum*delRatio;
	lg = 0;
	s = 0;
	while(s<delSum)
	{
		s += hist[lg];
		lg ++;
	}

	s = 0;
	rg = 255; 
	while(s<delSum)
	{
		s += hist[rg];
		rg --;
	}

	//3: 根据搜索的结果来对图像进行拉伸
	scale = 1 / (float)(rg-lg);
	for(i=0; i<imgSize; i++)
	{
		if(image[i]<=lg)
			image[i] = 0;
		else if(image[i]>=rg)
			image[i] = 255;
		else
			image[i] = (image[i]-lg)*scale*255;
	}
}



/* Geoimages Mosaic and color balance
*/
//void GeoMosaicWithBalance(vector<string> filesPath)
int GeoMosaicWithBalance(char** filesPath, int nFile, char* outpath, GDALDataType dataType)
{
	vector<stGeoInfo> geoArray; 

	GDALAllRegister();

	//if( filesPath.size()>1 )
	if(nFile>1)
	{
		int i,j;

		//get the area of mosaic
		//set the first image as the bench
		stGeoInfo geoinfo;
		GetGeoInformation(filesPath[0], geoinfo);
		geoArray.push_back(geoinfo);

		int zoneNumber = geoinfo.zoneNumber;
		//resolution
		double rx = geoinfo.dx;
		double ry = geoinfo.dy;
		double maxrx = rx;
		double maxry = fabs(ry);
		double minx,maxx,miny,maxy;
		minx = geoinfo.left;
		maxx = minx + geoinfo.wd*rx;
		maxy = geoinfo.top;
		miny = maxy + geoinfo.wd*ry;

		for(i=1; i<nFile; i++)
		{			
			GetGeoInformation(filesPath[i], geoinfo);
			geoArray.push_back(geoinfo);
			
			if(zoneNumber!=geoinfo.zoneNumber)
				continue;

			rx = geoinfo.dx;
			ry = geoinfo.dy;
			maxrx = max(rx, maxrx);
			maxry = max( fabs(ry), fabs(maxry) );
			printf("%lf %lf \n", rx, ry);

			minx = min(minx, geoinfo.left);
			maxx = max(maxx, geoinfo.left+geoinfo.wd*rx);
			miny = min(miny, geoinfo.top+geoinfo.ht*ry);
			maxy = max(maxy, geoinfo.top);
		}
		printf(" left:%lf right:%lf  top:%lf bottom:%lf \n", minx, maxx, maxy, miny);
	
		int mht = (maxy-miny) / fabs(geoArray[0].dy);
		int mwd = (maxx-minx) / geoArray[0].dx;
		printf("mosaic size: %d %d \n", mht, mwd);

		int nband = geoinfo.nband;
		
		//create mosaic image
		char **papszOptions = NULL;
		GDALDriver* poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
		GDALDataset* poDatasetMosaic = poDriver->Create(outpath, mwd, mht, nband, GDT_Byte, papszOptions );
		double  geoTransform[6];
		memset(geoTransform, 0, sizeof(double)*6);
		geoTransform[0] = minx;
		geoTransform[3] = maxy;
		geoTransform[1] = geoArray[0].dx;
		geoTransform[5] = geoArray[0].dy;
		poDatasetMosaic->SetGeoTransform(geoTransform);
		poDatasetMosaic->SetProjection(geoArray[0].projectRef);

		unsigned char* pBuffer = (unsigned char*)malloc( mht*mwd );
		memset(pBuffer, 0, mht*mwd);

		for(j=0; j<nband; j++)
		{
			printf("band %d \n", j);

			/*
			//the first image is set as the bench 
			GDALDataset* poDataset = (GDALDataset *) GDALOpen( filesPath[0], GA_ReadOnly );
			GDALRasterBand* poBand = poDataset->GetRasterBand( j+1 );
			int ht = geoArray[0].ht;
			int wd = geoArray[0].wd;
			unsigned char* pBench = (unsigned char*)malloc(ht*wd);
			poBand->RasterIO(GF_Read,0,0,wd,ht,pBench, wd, ht, GDT_Byte, 0, 0);
			GDALClose( (GDALDatasetH) poDataset );
			*/


			//get min and max value
			int minv = 10000000;
			int maxv = 0;
			for(i=0; i<nFile; i++)
			{
				//only deal with images within the same zone 
				if( zoneNumber!=geoArray[i].zoneNumber )
					continue;

				GDALDataset* poDataset = (GDALDataset *) GDALOpen( filesPath[i], GA_ReadOnly );
				GDALRasterBand* poBand = poDataset->GetRasterBand( j+1 );
				int ht = geoArray[i].ht;
				int wd = geoArray[i].wd;
				unsigned char* pImage = (unsigned char*)malloc(ht*wd);
				//int* pImage = (int*)malloc(ht*wd*sizeof(int));
				//assert(pImage!=NULL);
				poBand->RasterIO(GF_Read,0,0,wd,ht,pImage, wd, ht, dataType, 0, 0);
				GDALClose( (GDALDatasetH) poDataset );

				for(int k=0; k<ht*wd; k++)
				{
					if(minv>pImage[k])
						minv = pImage[k];
					if(maxv<pImage[k])
						maxv = pImage[k];
				}
			}

			//mosaic
			for(i=0; i<nFile; i++)
			{
				//only deal with images within the same zone 
				if( zoneNumber!=geoArray[i].zoneNumber )
					continue;
				
				GDALDataset* poDataset = (GDALDataset *) GDALOpen( filesPath[i], GA_ReadOnly );
				GDALRasterBand* poBand = poDataset->GetRasterBand( j+1 );
				int ht = geoArray[i].ht;
				int wd = geoArray[i].wd;
				unsigned char* pImage = (unsigned char*)malloc(ht*wd);
				//int* pImage = (int*)malloc(ht*wd*sizeof(int));
				//assert(pImage!=NULL);
				poBand->RasterIO(GF_Read,0,0,wd,ht,pImage, wd, ht, dataType, 0, 0);
				GDALClose( (GDALDatasetH) poDataset );
                
				/*
				//histogram matching 
				if(i!=0)
				{
					HistMatching(pBench, geoArray[0].ht, geoArray[0].wd, geoArray[0].wd, 
						pImage, ht, wd, wd);
				}
				*/

				//backward mapping
				int left = (geoArray[i].left - minx) / maxrx;
				int top  = (maxy - geoArray[i].top)  / maxry;
				int sw = geoArray[i].wd*geoArray[i].dx/maxrx;
				int sh = geoArray[i].ht*fabs(geoArray[i].dy)/maxry;
				printf("%d %d \n", left, top);

				for(int m=0; m<sh; m++)
					for(int n=0; n<sw; n++)
					{
						int row = (m*maxry)/fabs(geoArray[i].dy);
						int col = (n*maxrx)/geoArray[i].dx;
						col = max( 0, min(col,wd-1) );
						row = max( 0, min(row,ht-1) );
						
						if(pImage[row*wd+col]<3)
							continue;

						if( (top+m)>=mht || (left+n)>=mwd )
							continue;

						pBuffer[(top+m)*mwd+(left+n)] = (double)(pImage[row*wd+col]-minv) / (double)(maxv-minv)*255;
					}
				free(pImage);
			}

			//write data buffer
			//AutoLevelImage1(pBuffer, mht, mwd);
			GDALRasterBand* poBand = poDatasetMosaic->GetRasterBand( j+1 );
			poBand->RasterIO(GF_Write,0,0,mwd,mht,pBuffer,mwd,mht,GDT_Byte,0,0);
		}
		
		free(pBuffer);
		GDALClose( (GDALDatasetH) poDatasetMosaic );
	}
	
	//printf("Finished and the result is d:\\geomosaic.tif! \n");
	
	return 1;
}
