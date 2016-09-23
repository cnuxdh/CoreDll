#include "float.h"
#include "assert.h"

#include "geotiff.h"


//corelib
#include "CommonFuncs.h"


#include "gdal_priv.h"
#include "gdal.h"
#include "gdal_alg.h"
#include "ogr_spatialref.h"
#include "ogr_srs_api.h"
#include "cpl_string.h"
#include "cpl_conv.h"
#include "cpl_multiproc.h"
#include "cpl_port.h"


int ReadGeoFile(char* filename, stGeoInfo& geoData)
{
	FILE* fp = fopen(filename, "r");
	fscanf(fp, "%d ", &(geoData.zoneNumber));
	fscanf(fp, "%d %d ", &(geoData.ht), &(geoData.wd));
	fscanf(fp, "%lf %lf %lf ", &(geoData.left), &(geoData.top), &(geoData.dx));
	fclose(fp);

	geoData.dy = -geoData.dx;

	return 0;
}



int GetGeoInformation(const char* pszFilename, stGeoInfo& geoInfo)
{	
	//printf("[GetGeoInformation]  \n");
	
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");

	/*
	printf("GetDriver... \n");
	GDALDriver * poDriver = GetGDALDriverManager()->GetDriverByName("GTIFF");
	if( poDriver == NULL )
	{
		printf("Gdal Open file failed ! \n");
		return 0;
	}
	*/

	//printf("Open... %s \n", pszFilename);
	GDALDataset  *poDataset = NULL;
	poDataset = (GDALDataset *)GDALOpen( pszFilename, GA_ReadOnly );
	if( poDataset == NULL )
	{
		printf("Open failed! \n");
		return 0;
	}
	
	//read driver
	//printf( "Driver: %s/%s \n\n",
	//	poDataset->GetDriver()->GetDescription(), 
	//	poDataset->GetDriver()->GetMetadataItem( GDAL_DMD_LONGNAME ) );

	//read image size and channel
	//printf( "Size is %dx%dx%d \n\n", 
	//	poDataset->GetRasterXSize(), poDataset->GetRasterYSize(),
	//	poDataset->GetRasterCount() );

	//read projection
	int nLen = strlen( poDataset->GetProjectionRef() );
    geoInfo.projectRef = (const char*)malloc(nLen);
	strcpy( (char*)geoInfo.projectRef, poDataset->GetProjectionRef() );
	//geoInfo.projectRef = poDataset->GetProjectionRef();
	if( poDataset->GetProjectionRef()  != NULL )
	{
		//printf( "Projection is `%s'\n\n", poDataset->GetProjectionRef() );
	}
	

	//get the zone number
	OGRSpatialReference inpSR;
	char* Projection_str = (char*)(poDataset->GetProjectionRef());
	inpSR.importFromWkt(  &Projection_str );
	int zoneNumber =  inpSR.GetUTMZone();
	geoInfo.zoneNumber = zoneNumber;
	//printf("zone number: %d \n", zoneNumber);

	//read transform information
	double  adfGeoTransform[6];
	if( poDataset->GetGeoTransform( adfGeoTransform ) == CE_None )
	{
		//printf( "Origin = (%.6f,%.6f)\n",
		//	adfGeoTransform[0], adfGeoTransform[3] );
		//printf( "Pixel Size = (%.6f,%.6f)\n",
		//	adfGeoTransform[1], adfGeoTransform[5] );
	}

	//save geo-information into the structure
	geoInfo.left = adfGeoTransform[0];
	geoInfo.top  = adfGeoTransform[3];
	geoInfo.dx   = adfGeoTransform[1];
	geoInfo.dy   = adfGeoTransform[5];
	geoInfo.wd = poDataset->GetRasterXSize();
	geoInfo.ht = poDataset->GetRasterYSize();
	geoInfo.nband = poDataset->GetRasterCount();

	//Fetching a Raster Band
	GDALRasterBand  *poBand;
	int             nBlockXSize, nBlockYSize;
	int             bGotMin, bGotMax;
	double          adfMinMax[2];

	poBand = poDataset->GetRasterBand( 1 );
	poBand->GetBlockSize( &nBlockXSize, &nBlockYSize );
	
	//printf( "Block=%dx%d Type=%s, ColorInterp=%s \n\n",
	//	nBlockXSize, nBlockYSize,
	//	GDALGetDataTypeName(poBand->GetRasterDataType()),
	//	GDALGetColorInterpretationName(poBand->GetColorInterpretation()) );
		
	//poBand->GetBlockSize( &nBlockXSize, &nBlockYSize );
	//printf( "\n\n Block=%dx%d Type=%s, ColorInterp=%s \n\n",
	//	nBlockXSize, nBlockYSize,
	
	//printf("data type: %d \n", GDALGetDataTypeName( poBand->GetRasterDataType()) );
	//GDALGetColorInterpretationName(poBand->GetColorInterpretation()) );

	//adfMinMax[0] = poBand->GetMinimum( &bGotMin );
	//adfMinMax[1] = poBand->GetMaximum( &bGotMax );
	//if( ! (bGotMin && bGotMax) )
	//	GDALComputeRasterMinMax((GDALRasterBandH)poBand, TRUE, adfMinMax);
	//printf( "Min=%.3fd, Max=%.3f \n\n", adfMinMax[0], adfMinMax[1] );

	//if( poBand->GetOverviewCount() > 0 )
	//	printf( "Band has %d overviews.\n\n", poBand->GetOverviewCount() );
	//if( poBand->GetColorTable() != NULL )
	//	printf( "Band has a color table with %d entries.\n\n", 
	//	poBand->GetColorTable()->GetColorEntryCount() );


	//double minlon,maxlon,minlat,maxlat;
	//GroundToLL(poDataset, minlon, minlat, maxlon, maxlat);

	//close
	//printf("close dataset !\n");

	GDALClose( (GDALDatasetH)poDataset );

	//printf("reading over !\n");

	return 1;
}

void  GroundToLL(double gx, double gy, double& lon, double& lat, stGeoInfo geoInfo)
{
	/*
	double Geoform[20];
	int wd = geoInfo.wd; 
	int ht = geoInfo.ht; 
	double Geopos[4];
	Geopos[0] = geoInfo.left;   //top left X
	Geopos[1] = geoInfo.top;	//top left Y
	Geopos[2] = Geoform[0] + wd*geoInfo.dx;	// bottom right X
	Geopos[3] = Geoform[3] + ht*geoInfo.dy;	// bottom right Y	
	*/

	OGRSpatialReference inpSR;
	OGRSpatialReference otpSR;
	inpSR.importFromWkt((char **)(&geoInfo.projectRef));
	//printf("%s \n", Projection_str);
	otpSR.SetWellKnownGeogCS("WGS84");	
	OGRCoordinateTransformation *poTransform = OGRCreateCoordinateTransformation(&inpSR, &otpSR);	

	lon = gx;
	lat = gy;
	poTransform->Transform(1, &lon, &lat, NULL);//大地坐标到经纬度坐标

}



//ground coordinate to longitude/latitude using gdal
void GroundToLL(GDALDataset* poDataset, double& minlon, double& minlat, double& maxlon, double& maxlat)
{
	const char* Projection_str = poDataset->GetProjectionRef();
	
	double Geoform[20];
	poDataset->GetGeoTransform(Geoform);
	GDALRasterBand *poBand;
	poBand = poDataset->GetRasterBand(1);	
	int wd = poBand->GetXSize();
	int ht = poBand->GetYSize();
	double Geopos[4];
	Geopos[0] = Geoform[0];			// 左上角X
	Geopos[1] = Geoform[3];			// 左上角Y
	//Geopos[0] += fabs(Geoform[1])*BLOCK_SIZE;
	//Geopos[1] += fabs(Geoform[5])*BLOCK_SIZE;
	Geopos[2] = Geoform[0] + wd*Geoform[1] + ht*Geoform[2];	// 右下角X
	Geopos[3] = Geoform[3] + wd*Geoform[4] + ht*Geoform[5];	// 右下角Y	

	OGRSpatialReference inpSR;
	OGRSpatialReference otpSR;
	inpSR.importFromWkt((char **)(&Projection_str));
	printf("%s \n", Projection_str);
	otpSR.SetWellKnownGeogCS("WGS84");	
	OGRCoordinateTransformation *poTransform = OGRCreateCoordinateTransformation(&inpSR, &otpSR);	

	double lat[4],lon[4];
	lon[0] = Geopos[0];  lat[0] = Geopos[1];//上左
	lon[1] = Geopos[0] + fabs(Geoform[1])*wd;
	lat[1] = Geopos[1];//下左
	lon[2] = Geopos[0] + fabs(Geoform[1])*wd;
	lat[2] = Geopos[1] - fabs(Geoform[5])*ht;//上右
	lon[3] = Geopos[0]; 
	lat[3] = Geopos[1] - fabs(Geoform[5])*ht;//下右
	poTransform->Transform(4, lon, lat, NULL);//大地坐标到经纬度坐标

	//lontitude and latitude range
	minlon = 360;
	minlat = 180;
	maxlon = 0;
	maxlat = 0;
	for(int i=0; i<4; i++)
	{
		if(minlon>lon[i]) minlon=lon[i];
		if(maxlon<lon[i]) maxlon=lon[i];
		if(minlat>lat[i]) minlat=lat[i];
		if(maxlat<lat[i]) maxlat=lat[i];
	}
	printf(" lon: %lf %lf  lat: %lf %lf \n", minlon, maxlon, minlat, maxlat);
}

int PrintImageInfo(char* pszFilename)
{
	GDALDataset  *poDataset;
	GDALAllRegister();

	poDataset = (GDALDataset *) GDALOpen( pszFilename, GA_ReadOnly );
	if( poDataset == NULL )
	{
		printf("Open failed! \n");
		return 0;
	}

	//read driver
	printf( "Driver: %s/%s \n\n",
		poDataset->GetDriver()->GetDescription(), 
		poDataset->GetDriver()->GetMetadataItem( GDAL_DMD_LONGNAME ) );

	//read image size and channel
	printf( "Size is %dx%dx%d \n\n", 
		poDataset->GetRasterXSize(), poDataset->GetRasterYSize(),
		poDataset->GetRasterCount() );

	//read projection
	if( poDataset->GetProjectionRef()  != NULL )
		printf( "Projection is `%s'\n\n", poDataset->GetProjectionRef() );

	//read transform information
	double  adfGeoTransform[6];
	if( poDataset->GetGeoTransform( adfGeoTransform ) == CE_None )
	{
		printf( "Origin = (%.6f,%.6f)\n",
			adfGeoTransform[0], adfGeoTransform[3] );

		printf( "Pixel Size = (%.6f,%.6f)\n",
			adfGeoTransform[1], adfGeoTransform[5] );
	}



	//Fetching a Raster Band
	GDALRasterBand  *poBand;
	int             nBlockXSize, nBlockYSize;
	int             bGotMin, bGotMax;
	double          adfMinMax[2];

	poBand = poDataset->GetRasterBand( 1 );
	poBand->GetBlockSize( &nBlockXSize, &nBlockYSize );
	printf( "\n\n Block=%dx%d Type=%s, ColorInterp=%s \n\n",
		nBlockXSize, nBlockYSize,
		GDALGetDataTypeName(poBand->GetRasterDataType()),
		GDALGetColorInterpretationName(poBand->GetColorInterpretation()) );

	adfMinMax[0] = poBand->GetMinimum( &bGotMin );
	adfMinMax[1] = poBand->GetMaximum( &bGotMax );
	if( ! (bGotMin && bGotMax) )
		GDALComputeRasterMinMax((GDALRasterBandH)poBand, TRUE, adfMinMax);
	printf( "Min=%.3fd, Max=%.3f \n\n", adfMinMax[0], adfMinMax[1] );

	if( poBand->GetOverviewCount() > 0 )
		printf( "Band has %d overviews.\n\n", poBand->GetOverviewCount() );
	if( poBand->GetColorTable() != NULL )
		printf( "Band has a color table with %d entries.\n\n", 
		poBand->GetColorTable()->GetColorEntryCount() );


	double minlon,maxlon,minlat,maxlat;
	GroundToLL(poDataset, minlon, minlat, maxlon, maxlat);


	//close
	GDALClose( (GDALDatasetH) poDataset );

	return 1;
}

GDALDataType  GetDataType(const char* pszFilename)
{
	GDALDataset  *poDataset;

	GDALDataType ntype = GDT_Byte;

	GDALAllRegister();

	poDataset = (GDALDataset *) GDALOpen( pszFilename, GA_ReadOnly );
	if( poDataset == NULL )
	{
		printf("Open failed! \n");
		return ntype;
	}

	GDALRasterBand  *poBand;
	poBand = poDataset->GetRasterBand( 1 );
	ntype = poBand->GetRasterDataType();

	GDALClose( (GDALDatasetH) poDataset );	

	return ntype;
}



int  GdalWriteImageByteColor(char* savePath, unsigned char* r, unsigned char* g, unsigned char* b, int ht, int wd, stGeoInfo geoInfo)
{
	GDALAllRegister();
	GDALDriver* poDriver = NULL;
	GDALDataset *poDataset = NULL;   //GDAL数据集
	GDALRasterBand *poBand = NULL;
	char **papszOptions = NULL;
	//GDALRasterBand *poBand;
	poDriver = GetGDALDriverManager()->GetDriverByName("JPEG");
	if( poDriver == NULL )
	{
		printf("Gdal Open file failed ! \n");
		return 0;
	}
	poDataset = poDriver->Create(savePath, wd, ht, 3, GDT_Byte, papszOptions );
    
	double adfGeoTransform[6];
	adfGeoTransform[0] = geoInfo.left;
	adfGeoTransform[3] = geoInfo.top;
	adfGeoTransform[1] = geoInfo.dx;
	adfGeoTransform[5] = geoInfo.dy;
	poDataset->SetGeoTransform(adfGeoTransform);

	poBand = poDataset->GetRasterBand(1);
	poBand->RasterIO(GF_Write, 0, 0, wd, ht, r, wd, ht, GDT_Byte, 0, 0);
	poBand = poDataset->GetRasterBand(2);
	poBand->RasterIO(GF_Write, 0, 0, wd, ht, g, wd, ht, GDT_Byte, 0, 0);
	poBand = poDataset->GetRasterBand(3);
	poBand->RasterIO(GF_Write, 0, 0, wd, ht, b, wd, ht, GDT_Byte, 0, 0);

	//close
	GDALClose( (GDALDatasetH) poDataset );
	return 1;


}

int  GdalWriteImageByteColor(char* savePath, unsigned char* r, unsigned char* g, unsigned char* b, int ht, int wd)
{
	GDALAllRegister();

	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");

	char **papszOptions = NULL;

	char* drvName="MEM";
	GDALDriver *memDriver;
	memDriver = GetGDALDriverManager()->GetDriverByName(drvName);
	GDALDataset* ds_preview = memDriver->Create("default", wd, ht, 3, GDT_Byte, papszOptions);

	GDALRasterBand* ba_preview=ds_preview->GetRasterBand(1);
	ba_preview->RasterIO(GF_Write, 0, 0, wd, ht, r, wd, ht, GDT_Byte, 0, 0);
	ba_preview=ds_preview->GetRasterBand(2);
	ba_preview->RasterIO(GF_Write, 0, 0, wd, ht, g, wd, ht, GDT_Byte, 0, 0);
	ba_preview=ds_preview->GetRasterBand(3);
	ba_preview->RasterIO(GF_Write, 0, 0, wd, ht, b, wd, ht, GDT_Byte, 0, 0);
	//ba_preview=ds_preview->GetRasterBand(4);
	//ba_preview->RasterIO(GF_Write, 0, 0, wd, ht, b, wd, ht, GDT_Byte, 0, 0);

	char* jpgDrvName="JPEG";
	GDALDriver *jpgDriver=GetGDALDriverManager()->GetDriverByName(jpgDrvName);
	GDALDataset *dsjpg=jpgDriver->CreateCopy(savePath, ds_preview,0,NULL,NULL,NULL);

	//close
	GDALClose( (GDALDatasetH) ds_preview );
	GDALClose( (GDALDatasetH) dsjpg );
	
	return 1;
}

int  GdalWriteColorImageUShort(char* savePath, 
	unsigned short* r, unsigned short* g, unsigned short* b, int ht, int wd,
	stGeoInfo geoInfo)
{

	GDALAllRegister();
	GDALDriver* poDriver = NULL;
	GDALDataset *poDataset = NULL;   //GDAL数据集
	GDALRasterBand *poBand = NULL;
	char **papszOptions = NULL;

	//GDALRasterBand *poBand;
	poDriver = GetGDALDriverManager()->GetDriverByName("GTIFF");
	if( poDriver == NULL )
	{
		printf("Gdal Open file failed ! \n");
		return 0;
	}

	poDataset = poDriver->Create(savePath, wd, ht, 3, GDT_UInt16, papszOptions );

	double  geoTransform[6];
	memset(geoTransform, 0, sizeof(double)*6);
	geoTransform[0] = geoInfo.left;
	geoTransform[3] = geoInfo.top;
	geoTransform[1] = geoInfo.dx;
	geoTransform[5] = geoInfo.dy;
	poDataset->SetGeoTransform(geoTransform);
	poDataset->SetProjection(geoInfo.projectRef);

	poBand = poDataset->GetRasterBand(1);
	poBand->RasterIO(GF_Write, 0, 0, wd, ht, r, wd, ht, GDT_UInt16, 0, 0);
	poBand = poDataset->GetRasterBand(2);
	poBand->RasterIO(GF_Write, 0, 0, wd, ht, g, wd, ht, GDT_UInt16, 0, 0);
	poBand = poDataset->GetRasterBand(3);
	poBand->RasterIO(GF_Write, 0, 0, wd, ht, b, wd, ht, GDT_UInt16, 0, 0);

	//close
	GDALClose( (GDALDatasetH) poDataset );

	return 0;
}



int  GdalWriteImageColor(char* savePath, unsigned char* r, unsigned char* g, unsigned char* b, 
	int ht, int wd, stGeoInfo geoInfo)
{
	GDALAllRegister();
	GDALDriver* poDriver = NULL;
	GDALDataset *poDataset = NULL;   //GDAL数据集
	GDALRasterBand *poBand = NULL;
	char **papszOptions = NULL;
	//GDALRasterBand *poBand;
	poDriver = GetGDALDriverManager()->GetDriverByName("GTIFF");
	if( poDriver == NULL )
	{
		printf("Gdal Open file failed ! \n");
		return 0;
	}

	poDataset = poDriver->Create(savePath, wd, ht, 3, GDT_Byte, papszOptions );

	double  geoTransform[6];
	memset(geoTransform, 0, sizeof(double)*6);
	geoTransform[0] = geoInfo.left;
	geoTransform[3] = geoInfo.top;
	geoTransform[1] = geoInfo.dx;
	geoTransform[5] = geoInfo.dy;
	poDataset->SetGeoTransform(geoTransform);
	poDataset->SetProjection(geoInfo.projectRef);
	
	poBand = poDataset->GetRasterBand(1);
	poBand->RasterIO(GF_Write, 0, 0, wd, ht, r, wd, ht, GDT_Byte, 0, 0);
	poBand = poDataset->GetRasterBand(2);
	poBand->RasterIO(GF_Write, 0, 0, wd, ht, g, wd, ht, GDT_Byte, 0, 0);
	poBand = poDataset->GetRasterBand(3);
	poBand->RasterIO(GF_Write, 0, 0, wd, ht, b, wd, ht, GDT_Byte, 0, 0);

	//close
	GDALClose( (GDALDatasetH) poDataset );
	return 1;
}

int  GdalWriteByte(char* savePath, unsigned char* pBuffer, int ht, int wd, stGeoInfo geoInfo)
{

	GDALAllRegister();
	GDALDriver* poDriver   = NULL;
	GDALDataset *poDataset = NULL;   //GDAL数据集
	GDALRasterBand *poBand = NULL;
	char **papszOptions    = NULL;
	//GDALRasterBand *poBand;
	poDriver = GetGDALDriverManager()->GetDriverByName("GTIFF");
	if( poDriver == NULL )
	{
		printf("Gdal Open file failed ! \n");
		return 0;
	}

	poDataset = poDriver->Create(savePath, wd, ht, 1, GDT_Byte, papszOptions );

	double  geoTransform[6];
	memset(geoTransform, 0, sizeof(double)*6);
	geoTransform[0] = geoInfo.left;
	geoTransform[3] = geoInfo.top;
	geoTransform[1] = geoInfo.dx;
	geoTransform[5] = geoInfo.dy;
	poDataset->SetGeoTransform(geoTransform);
	poDataset->SetProjection(geoInfo.projectRef);

	poBand = poDataset->GetRasterBand(1);
	poBand->RasterIO(GF_Write, 0, 0, wd, ht, pBuffer, wd, ht, GDT_Byte, 0, 0);

	//close
	GDALClose( (GDALDatasetH) poDataset );

	return 1;

}


int  GdalWriteUShort(char* savePath, unsigned short* pBuffer, int ht, int wd, stGeoInfo geoInfo)
{

	GDALAllRegister();
	GDALDriver* poDriver   = NULL;
	GDALDataset *poDataset = NULL;   //GDAL数据集
	GDALRasterBand *poBand = NULL;
	char **papszOptions    = NULL;
	
	//GDALRasterBand *poBand;
	poDriver = GetGDALDriverManager()->GetDriverByName("GTIFF");
	if( poDriver == NULL )
	{
		printf("Gdal Open file failed ! \n");
		return 0;
	}

	poDataset = poDriver->Create(savePath, wd, ht, 1, GDT_UInt16, papszOptions );

	double  geoTransform[6];
	memset(geoTransform, 0, sizeof(double)*6);
	geoTransform[0] = geoInfo.left;
	geoTransform[3] = geoInfo.top;
	geoTransform[1] = geoInfo.dx;
	geoTransform[5] = geoInfo.dy;
	poDataset->SetGeoTransform(geoTransform);
	poDataset->SetProjection(geoInfo.projectRef);

	poBand = poDataset->GetRasterBand(1);
	poBand->RasterIO(GF_Write, 0, 0, wd, ht, pBuffer, wd, ht, GDT_UInt16, 0, 0);

	//close
	GDALClose( (GDALDatasetH) poDataset );

	return 1;

}

int  GdalWriteUInt(char* savePath, unsigned int* pBuffer, int ht, int wd, stGeoInfo geoInfo)
{
	GDALAllRegister();
	GDALDriver* poDriver   = NULL;
	GDALDataset *poDataset = NULL;   //GDAL数据集
	GDALRasterBand *poBand = NULL;
	char **papszOptions    = NULL;
	//GDALRasterBand *poBand;
	poDriver = GetGDALDriverManager()->GetDriverByName("GTIFF");
	if( poDriver == NULL )
	{
		printf("Gdal Open file failed ! \n");
		return 0;
	}

	poDataset = poDriver->Create(savePath, wd, ht, 1, GDT_UInt32, papszOptions );

	double  geoTransform[6];
	memset(geoTransform, 0, sizeof(double)*6);
	geoTransform[0] = geoInfo.left;
	geoTransform[3] = geoInfo.top;
	geoTransform[1] = geoInfo.dx;
	geoTransform[5] = geoInfo.dy;
	poDataset->SetGeoTransform(geoTransform);
	poDataset->SetProjection(geoInfo.projectRef);

	poBand = poDataset->GetRasterBand(1);
	poBand->RasterIO(GF_Write, 0, 0, wd, ht, pBuffer, wd, ht, GDT_UInt32, 0, 0);

	//close
	GDALClose( (GDALDatasetH) poDataset );

	return 1;

}

int  GdalWriteInt(char* savePath, int* pBuffer, int ht, int wd, stGeoInfo geoInfo)
{
	GDALAllRegister();
	GDALDriver* poDriver   = NULL;
	GDALDataset *poDataset = NULL;   //GDAL数据集
	GDALRasterBand *poBand = NULL;
	char **papszOptions    = NULL;
	//GDALRasterBand *poBand;
	poDriver = GetGDALDriverManager()->GetDriverByName("GTIFF");
	if( poDriver == NULL )
	{
		printf("Gdal Open file failed ! \n");
		return 0;
	}

	poDataset = poDriver->Create(savePath, wd, ht, 1, GDT_Int32, papszOptions );

	double  geoTransform[6];
	memset(geoTransform, 0, sizeof(double)*6);
	geoTransform[0] = geoInfo.left;
	geoTransform[3] = geoInfo.top;
	geoTransform[1] = geoInfo.dx;
	geoTransform[5] = geoInfo.dy;
	poDataset->SetGeoTransform(geoTransform);
	poDataset->SetProjection(geoInfo.projectRef);

	poBand = poDataset->GetRasterBand(1);
	poBand->RasterIO(GF_Write, 0, 0, wd, ht, pBuffer, wd, ht, GDT_Int32, 0, 0);

	//close
	GDALClose( (GDALDatasetH) poDataset );

	return 1;

}

int  GdalWriteFloat(char* savePath, float* pBuffer, int ht, int wd, stGeoInfo geoInfo)
{
	GDALAllRegister();
	GDALDriver* poDriver   = NULL;
	GDALDataset *poDataset = NULL;   //GDAL数据集
	GDALRasterBand *poBand = NULL;
	char **papszOptions    = NULL;
	//GDALRasterBand *poBand;
	poDriver = GetGDALDriverManager()->GetDriverByName("GTIFF");
	if( poDriver == NULL )
	{
		printf("Gdal Open file failed ! \n");
		return 0;
	}

	poDataset = poDriver->Create(savePath, wd, ht, 1, GDT_Float32, papszOptions );

	double  geoTransform[6];
	memset(geoTransform, 0, sizeof(double)*6);
	geoTransform[0] = geoInfo.left;
	geoTransform[3] = geoInfo.top;
	geoTransform[1] = geoInfo.dx;
	geoTransform[5] = geoInfo.dy;
	poDataset->SetGeoTransform(geoTransform);
	
	if(geoInfo.projectRef != NULL)
		poDataset->SetProjection(geoInfo.projectRef);

	poBand = poDataset->GetRasterBand(1);
	poBand->RasterIO(GF_Write, 0, 0, wd, ht, pBuffer, wd, ht, GDT_Float32, 0, 0);

	//close
	GDALClose( (GDALDatasetH) poDataset );

	return 1;
}


//write buffer into image file using Gdal lib
int GdalWriteImageByte(char* savePath, unsigned char* pBuffer, int ht, int wd)
{	
	GDALAllRegister();
	GDALDriver* poDriver = NULL;
	GDALDataset *poDataset = NULL;   //GDAL数据集
	GDALRasterBand *poBand = NULL;
	char **papszOptions = NULL;
	//GDALRasterBand *poBand;
	poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
	if( poDriver == NULL )
	{
		printf("Gdal Open file failed ! \n");
		return 0;
	}
	poDataset = poDriver->Create(savePath, wd, ht, 1, GDT_Byte, papszOptions );
	poBand = poDataset->GetRasterBand(1);
	poBand->RasterIO(GF_Write, 0, 0, wd, ht, pBuffer, wd, ht, GDT_Byte, 0, 0);
	//close
	GDALClose( (GDALDatasetH) poDataset );
	return 1;
}


//write buffer into image file using Gdal lib
int GdalWriteImageUShort(char* savePath, unsigned short* pBuffer, int ht, int wd)
{	
	GDALAllRegister();
	GDALDriver* poDriver = NULL;
	GDALDataset *poDataset = NULL;   //GDAL数据集
	GDALRasterBand *poBand = NULL;
	char **papszOptions = NULL;
	//GDALRasterBand *poBand;
	
	poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
	if( poDriver == NULL )
	{
		printf("Gdal Open file failed ! \n");
		return 0;
	}
	poDataset = poDriver->Create(savePath, wd, ht, 1, GDT_UInt16, papszOptions );
	poBand = poDataset->GetRasterBand(1);
	poBand->RasterIO(GF_Write, 0, 0, wd, ht, pBuffer, wd, ht, GDT_UInt16, 0, 0);
	//close
	GDALClose( (GDALDatasetH) poDataset );
	return 1;
}

//jpg is different from tif, and must be created by "CreateCopy"
int GdalWriteJpgCopy(char* savePath, unsigned char* pBuffer, int ht, int wd)
{	
	GDALAllRegister();
	char **papszOptions = NULL;

	char* drvName="MEM";
	GDALDriver *memDriver;
	memDriver = GetGDALDriverManager()->GetDriverByName(drvName);
	GDALDataset* ds_preview = memDriver->Create("default", wd, ht, 1, GDT_Byte, papszOptions);

	GDALRasterBand* ba_preview=ds_preview->GetRasterBand(1);
	ba_preview->RasterIO(GF_Write, 0, 0, wd, ht, pBuffer, wd, ht, GDT_Byte, 0, 0);

	char* jpgDrvName="JPEG";
	GDALDriver *jpgDriver=GetGDALDriverManager()->GetDriverByName(jpgDrvName);
	GDALDataset *dsjpg=jpgDriver->CreateCopy(savePath, ds_preview,0,NULL,NULL,NULL);

	//close
	GDALClose( (GDALDatasetH) ds_preview );
	GDALClose( (GDALDatasetH) dsjpg );

	return 1;
}


void GetImageProject(GDALDataset* poDataset, double* Geopos )
{	
	const char* Projection_str = poDataset->GetProjectionRef();
	double Geoform[20];
	poDataset->GetGeoTransform(Geoform);
	GDALRasterBand *poBand;
	poBand = poDataset->GetRasterBand(1);	
	int wd = poBand->GetXSize();
	int ht = poBand->GetYSize();
	//double Geopos[4];
	Geopos[0] = Geoform[0];			// 左上角X
	Geopos[1] = Geoform[3];			// 左上角Y
	//Geopos[0] += fabs(Geoform[1])*BLOCK_SIZE;
	//Geopos[1] += fabs(Geoform[5])*BLOCK_SIZE;
	Geopos[2] = Geoform[0] + wd*Geoform[1] + ht*Geoform[2];	// 右下角X
	Geopos[3] = Geoform[3] + wd*Geoform[4] + ht*Geoform[5];	// 右下角Y		
}
/*
//ground coordinate to longitude/latitude using gdal
void GroundToLL(GDALDataset* poDataset, double& minlon, double& minlat, double& maxlon, double& maxlat)
{
	const char* Projection_str = poDataset->GetProjectionRef();
	double Geoform[20];
	poDataset->GetGeoTransform(Geoform);
	GDALRasterBand *poBand;
	poBand = poDataset->GetRasterBand(1);	
	int wd = poBand->GetXSize();
	int ht = poBand->GetYSize();
	double Geopos[4];
	Geopos[0] = Geoform[0];			// 左上角X
	Geopos[1] = Geoform[3];			// 左上角Y
	//Geopos[0] += fabs(Geoform[1])*BLOCK_SIZE;
	//Geopos[1] += fabs(Geoform[5])*BLOCK_SIZE;
	Geopos[2] = Geoform[0] + wd*Geoform[1] + ht*Geoform[2];	// 右下角X
	Geopos[3] = Geoform[3] + wd*Geoform[4] + ht*Geoform[5];	// 右下角Y	

	OGRSpatialReference inpSR;
	OGRSpatialReference otpSR;
	inpSR.importFromWkt((char **)(&Projection_str));
	printf("%s \n", Projection_str);
	otpSR.SetWellKnownGeogCS("WGS84");	
	OGRCoordinateTransformation *poTransform = OGRCreateCoordinateTransformation(&inpSR, &otpSR);	
		
	double lat[4],lon[4];
	lon[0] = Geopos[0];  lat[0] = Geopos[1];//上左
	lon[1] = Geopos[0] + fabs(Geoform[1])*wd;
	lat[1] = Geopos[1];//下左
	lon[2] = Geopos[0] + fabs(Geoform[1])*wd;
	lat[2] = Geopos[1] - fabs(Geoform[5])*ht;//上右
	lon[3] = Geopos[0]; 
	lat[3] = Geopos[1] - fabs(Geoform[5])*ht;//下右
	poTransform->Transform(4, lon, lat, NULL);//大地坐标到经纬度坐标

	//lontitude and latitude range
	minlon = 360;
	minlat = 180;
	maxlon = 0;
	maxlat = 0;
	for(int i=0; i<4; i++)
	{
		if(minlon>lon[i]) minlon=lon[i];
		if(maxlon<lon[i]) maxlon=lon[i];
		if(minlat>lat[i]) minlat=lat[i];
		if(maxlat<lat[i]) maxlat=lat[i];
	}
}
*/


/* 将数据写入文件，投影方式采用经纬度方式   
*/
int  GdalWriteFloatLL(char* savePath, float* pBuffer, int ht, int wd, 
					  double minlon, double maxlon, double minlax, double maxlax)
{
	char *projection = "GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563,AUTHORITY[\"EPSG\",\"7030\"]],AUTHORITY[\"EPSG\",\"6326\"]],PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433],AUTHORITY[\"EPSG\",\"4326\"]]";
	double latlonPoint[6];

	double pixelsize = (maxlon-minlon) / (double)(wd);

	latlonPoint[0] = minlon; 
	latlonPoint[1] = pixelsize;
	latlonPoint[2] = 0;
	latlonPoint[3] = maxlax; 
	latlonPoint[4] = 0;
	latlonPoint[5] = -pixelsize;	

	char* tifDrvName="GTiff";
	char **optionstif=NULL;

	GDALAllRegister();

	GDALDriver *tifDriver=GetGDALDriverManager()->GetDriverByName(tifDrvName);
	GDALDataset *dstif;
	GDALRasterBand *baTif;

	dstif = tifDriver->Create(savePath,wd,ht,1,GDT_Float32,optionstif);
	baTif = dstif->GetRasterBand(1);
	baTif->RasterIO(GF_Write, 0, 0, wd, ht, 
		pBuffer, wd, ht, GDT_Float32, 0, 0);

	dstif->SetProjection(projection);
	dstif->SetGeoTransform(latlonPoint);

	GDALClose((GDALDatasetH)dstif);
	return 1;
}

int  GdalWriteFloatLL(char* savePath, 
								 unsigned char* r, unsigned char* g, unsigned char* b,
								 int ht, int wd, 
								 double minlon, double maxlon, double minlax, double maxlax)
{

	char *projection = "GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563,AUTHORITY[\"EPSG\",\"7030\"]],AUTHORITY[\"EPSG\",\"6326\"]],PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433],AUTHORITY[\"EPSG\",\"4326\"]]";
	double latlonPoint[6];

	double pixelsize = (maxlon-minlon) / (double)(wd);

	latlonPoint[0] = minlon; 
	latlonPoint[1] = pixelsize;
	latlonPoint[2] = 0;
	latlonPoint[3] = maxlax; 
	latlonPoint[4] = 0;
	latlonPoint[5] = -pixelsize;	

	char* tifDrvName="GTiff";
	char **optionstif=NULL;

	GDALAllRegister();

	GDALDriver *tifDriver=GetGDALDriverManager()->GetDriverByName(tifDrvName);
	GDALDataset *dstif;
	GDALRasterBand *baTif;

	dstif = tifDriver->Create(savePath,wd,ht,3,GDT_Byte,optionstif);
	baTif = dstif->GetRasterBand(1);
	baTif->RasterIO(GF_Write, 0, 0, wd, ht, 
		r, wd, ht, GDT_Byte, 0, 0);

	baTif = dstif->GetRasterBand(2);
	baTif->RasterIO(GF_Write, 0, 0, wd, ht, 
		g, wd, ht, GDT_Byte, 0, 0);

	baTif = dstif->GetRasterBand(3);
	baTif->RasterIO(GF_Write, 0, 0, wd, ht, 
		b, wd, ht, GDT_Byte, 0, 0);

	dstif->SetProjection(projection);
	dstif->SetGeoTransform(latlonPoint);

	GDALClose((GDALDatasetH)dstif);
	return 1;
}


/* Mosaic the files including Geoinfo such as GeoTiff,IMG, etc.
   written by xiedonghai, 2013.8.13
*/
void GeoImageMosaic(char** files, int nfile)
{
	/*
	//get the area of each file
	for(int i=0; i<nfile; i++)
	{
		stGeoInfo gi;
		printf("%s \n", files[i]);
		GetGeoInformation(files[i], gi);
	}
	*/

	GDALDataset  *poDataset;
	GDALAllRegister();

	poDataset = (GDALDataset *) GDALOpen( files[0], GA_ReadOnly );
	if( poDataset == NULL )
	{
		printf("Open failed! \n");
		return ;
	}
}


//write buffer into image file using Gdal lib
int GdalWriteImageInt(char* savePath, short* pBuffer, int ht, int wd, 
					   const char* projectRef, double* geoTransform)
{	
	GDALAllRegister();
	GDALDriver* poDriver = NULL;
	GDALDataset *poDataset = NULL;   //GDAL数据集
	GDALRasterBand *poBand = NULL;
	char **papszOptions = NULL;
	//GDALRasterBand *poBand;

	poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
	if( poDriver == NULL )
	{
		printf("Gdal Open file failed ! \n");
		return 0;
	}
	poDataset = poDriver->Create(savePath, wd, ht, 1, GDT_Int16, papszOptions );
	poDataset->SetGeoTransform(geoTransform);
	poDataset->SetProjection(projectRef);
	poBand = poDataset->GetRasterBand(1);
	poBand->RasterIO(GF_Write, 0, 0, wd, ht, pBuffer, wd, ht, GDT_Int16, 0, 0);
	//close
	GDALClose( (GDALDatasetH) poDataset );
	return 1;
}

int GeotiffMedianFilter(char* filename)
{
	GDALAllRegister();

	GDALDataset  *poDataset = (GDALDataset *) GDALOpen( filename, GA_ReadOnly );
	if( poDataset == NULL )
	{
		printf("Open failed! \n");
		return 0;
	}
	
	//read image buffer
	GDALRasterBand *poBand = poDataset->GetRasterBand(1);
	GDALDataType nType = poBand->GetRasterDataType();
	int wd = poBand->GetXSize();
	int ht = poBand->GetYSize();
	short* pBuffer = (short*)malloc(ht*wd*sizeof(short));		
	poBand->RasterIO(GF_Read,0,0,wd,ht,pBuffer,wd,ht,GDT_Int16,0,0);

	//median filter
	short* pRes = (short*)malloc(ht*wd*sizeof(short));		
	memcpy(pRes, pBuffer, sizeof(short)*ht*wd);

	short neibor[9];
	int i,j,m,n;
	for( j=1; j<ht-1; j++)
		for( i=1; i<wd-1; i++)
		{
			if(pBuffer[j*wd+i]==0)
			{
				int index = 0;
				for( m=-1; m<=1; m++)
					for( n=-1; n<=1; n++)
					{
						neibor[index] = pBuffer[ (j+m)*wd + (i+n) ];
						index++;
					}
				//sort
				for(m=0; m<9; m++)
					for(n=m+1; n<9; n++)
					{
						if(neibor[m]<neibor[n])
						{
							int t = neibor[m];
							neibor[m] = neibor[n];
							neibor[n] = t;
						}
					}
				//replace
				pRes[j*wd+i] = neibor[4];
			}
		}
	//memcpy(pBuffer, pRes, sizeof(short)*ht*wd);


	for( j=0; j<ht*wd; j++)
	{
		if(pRes[j]<0)
			pRes[j] = 0;
	}

	//overwrite the file
	double  adfGeoTransform[6];
	poDataset->GetGeoTransform( adfGeoTransform );
	const char* projectRef = poDataset->GetProjectionRef();
	GdalWriteImageInt( "d:\\dem.tif", pRes, ht, wd, projectRef, adfGeoTransform );	
	
	free(pBuffer);
	free(pRes);
	GDALClose( (GDALDatasetH) poDataset );
	return 1;
}


int  ReadGeoFileShort(char* filePath, int bandId, short** pBuffer, int& ht, int& wd)
{
	GDALAllRegister();

	GDALDataset *poDataset = (GDALDataset*)GDALOpen(filePath, GA_ReadOnly); 
	GDALRasterBand *poBand  = poDataset->GetRasterBand(bandId);
	wd = poBand->GetXSize();
	ht = poBand->GetYSize();
	//*pBuffer=(unsigned char*)CPLMalloc(sizeof(unsigned char)*wd*ht);
	*pBuffer=(short*)malloc(sizeof(short)*wd*ht);
	poBand->RasterIO(GF_Read,0,0,wd,ht,*pBuffer,wd,ht,GDT_Int16,0,0);

	GDALClose( (GDALDatasetH) poDataset );

	return 1;
}



int  ReadGeoFileUShort(char* filePath, int bandId, unsigned short** pBuffer, int& ht, int& wd)
{
	GDALAllRegister();

	GDALDataset *poDataset = (GDALDataset*)GDALOpen(filePath, GA_ReadOnly); 
	GDALRasterBand *poBand  = poDataset->GetRasterBand(bandId);
	wd = poBand->GetXSize();
	ht = poBand->GetYSize();
	//*pBuffer=(unsigned char*)CPLMalloc(sizeof(unsigned char)*wd*ht);
	*pBuffer=(unsigned short*)malloc(sizeof(unsigned short)*wd*ht);
	poBand->RasterIO(GF_Read,0,0,wd,ht,*pBuffer,wd,ht,GDT_UInt16,0,0);

	GDALClose( (GDALDatasetH) poDataset );

	return 1;
}

int  ReadGeoFileFloat(char* filePath, int bandId, float** pBuffer, int& ht, int& wd)
{
	GDALAllRegister();

	GDALDataset *poDataset = (GDALDataset*)GDALOpen(filePath, GA_ReadOnly); 
	GDALRasterBand *poBand  = poDataset->GetRasterBand(bandId);
	wd = poBand->GetXSize();
	ht = poBand->GetYSize();

	*pBuffer=(float*)malloc(sizeof(float)*wd*ht);
	if(pBuffer==NULL)
	{
		printf("Memory malloc failed! \n");
	}

	poBand->RasterIO(GF_Read,0,0,wd,ht,*pBuffer,wd,ht,GDT_Float32,0,0);

	GDALClose( (GDALDatasetH) poDataset );

	return 1;
}



int  ReadGeoFileInt(char* filePath, int bandId, int** pBuffer, int& ht, int& wd)
{
	GDALAllRegister();

	GDALDataset *poDataset = (GDALDataset*)GDALOpen(filePath, GA_ReadOnly); 
	GDALRasterBand *poBand  = poDataset->GetRasterBand(bandId);
	wd = poBand->GetXSize();
	ht = poBand->GetYSize();

	*pBuffer = (int*)malloc(sizeof(int)*wd*ht);
	poBand->RasterIO(GF_Read,0,0,wd,ht,*pBuffer,wd,ht,GDT_Int32,0,0);

	GDALClose( (GDALDatasetH) poDataset );

	return 1;
}

int  ReadGeoFileUInt(char* filePath, int bandId, unsigned int** pBuffer, int& ht, int& wd)
{
	GDALAllRegister();

	GDALDataset *poDataset = (GDALDataset*)GDALOpen(filePath, GA_ReadOnly); 
	GDALRasterBand *poBand  = poDataset->GetRasterBand(bandId);
	wd = poBand->GetXSize();
	ht = poBand->GetYSize();

	*pBuffer = (unsigned int*)malloc(sizeof(unsigned int)*wd*ht);
	poBand->RasterIO(GF_Read,0,0,wd,ht,*pBuffer,wd,ht,GDT_UInt32,0,0);

	GDALClose( (GDALDatasetH) poDataset );

	return 1;
}


int ReadGeoFileGeneral(char* filepath, int bandId, void* pBuffer, GDALDataType nType)
{
	GDALAllRegister();

	GDALDataset *poDataset = (GDALDataset*)GDALOpen(filepath, GA_ReadOnly); 
	GDALRasterBand *poBand  = poDataset->GetRasterBand(bandId);
	
	int wd = poBand->GetXSize();
	int ht = poBand->GetYSize();

	poBand->RasterIO(GF_Read, 0, 0, wd, ht, pBuffer, wd, ht, nType,0,0);

	GDALClose( (GDALDatasetH) poDataset );

	return 0;
}

DLL_EXPORT int  ReadGeoFileByte(char* filePath, int bandId, double ratio, 
								unsigned char** pBuffer, int& ht, int& wd)
{
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");

	GDALDataset *poDataset = (GDALDataset*)GDALOpen(filePath, GA_ReadOnly); 
	GDALRasterBand *poBand  = poDataset->GetRasterBand(bandId);
	int swd = poBand->GetXSize();
	int sht = poBand->GetYSize();

	ht = sht*ratio;
	wd = swd*ratio;

	int ntype = poBand->GetRasterDataType();
	*pBuffer=(unsigned char*)malloc(sizeof(unsigned char)*wd*ht);

	unsigned short* pUShort = NULL;
	short* pShort = NULL;
	unsigned int* pUInt = NULL;
	int* pInt = NULL; 
	float* pFloat = NULL;
	double* pDouble = NULL;
	double dMin, dMax;

	switch(ntype)
	{
	case 1: //byte
		poBand->RasterIO(GF_Read,0,0,swd,sht,*pBuffer,wd,ht,GDT_Byte,0,0);
		break;
	case 2: //16 bit unsigned int
		pUShort = (unsigned short*)malloc(ht*wd*2);
		poBand->RasterIO(GF_Read,0,0,swd,sht,pUShort,wd,ht,GDT_UInt16,0,0);
		GetArrayMinMax(pUShort, wd*ht, dMin, dMax);
		for(int i=0; i<wd*ht; i++)
		{
			if( pUShort[i]==0 )
				continue;
			(*pBuffer)[i] = (pUShort[i] - dMin) / (dMax-dMin)*255;
		}
		break;
	case 3: //16 bit signed int
		pShort = (short*)malloc(ht*wd*2);
		poBand->RasterIO(GF_Read,0,0,swd,sht,pShort,wd,ht,GDT_Int16,0,0);
		GetArrayMinMax(pShort, wd*ht, dMin, dMax);
		for(int i=0; i<wd*ht; i++)
		{
			if( pShort[i]==0 )
				continue;
			(*pBuffer)[i] = (pShort[i] - dMin) / (dMax-dMin)*255;
		}
		break;
	case 4: //32 bit unsigned int
		pUInt = (unsigned int*)malloc(ht*wd*4);
		poBand->RasterIO(GF_Read,0,0,swd,sht,pUInt,wd,ht,GDT_UInt32,0,0);
		GetArrayMinMax(pUInt, wd*ht, dMin, dMax);
		for(int i=0; i<wd*ht; i++)
		{
			if(pUInt[i]==0)
				continue;
			(*pBuffer)[i] = (pUInt[i] - dMin) / (dMax-dMin)*255;
		}
		break;
	case 5: //32 bit signed int
		pInt = (int*)malloc(ht*wd*4);
		poBand->RasterIO(GF_Read,0,0,swd,sht,pInt,wd,ht,GDT_Int32,0,0);
		GetArrayMinMax(pInt, wd*ht, dMin, dMax);
		for(int i=0; i<wd*ht; i++)
		{
			if( pInt[i]==0 )
				continue;
			(*pBuffer)[i] = (pInt[i] - dMin) / (dMax-dMin)*255;
		}
		break;
	case 6: //32 bit float
		pFloat = (float*)malloc(ht*wd*4);
		poBand->RasterIO(GF_Read,0,0,swd,sht,pUInt,wd,ht,GDT_Float32,0,0);
		GetArrayMinMax(pFloat, wd*ht, dMin, dMax);
		for(int i=0; i<wd*ht; i++)
		{
			if(pFloat[i]==0)
				continue;
			(*pBuffer)[i] = (pFloat[i] - dMin) / (dMax-dMin)*255;
		}
		break;
	case 7: //64 bit float
		pDouble = (double*)malloc(ht*wd*8);
		poBand->RasterIO(GF_Read,0,0,swd,sht,pDouble,wd,ht,GDT_Float64,0,0);
		GetArrayMinMax(pDouble, wd*ht, dMin, dMax);
		for(int i=0; i<wd*ht; i++)
		{
			if(pDouble[i]==0)
				continue;
			(*pBuffer)[i] = (pDouble[i] - dMin) / (dMax-dMin)*255;
		}
		break;
	default:
		break;		
	}

	free(pUShort);
	free(pShort);
	free(pUInt);
	free(pInt);
	free(pFloat);
	free(pDouble);

	GDALClose( (GDALDatasetH) poDataset );

	return 1;
}


int  ReadGeoFileByte(char* filePath, int bandId, unsigned char** pBuffer, int& ht, int& wd)
{
	GDALAllRegister();

	GDALDataset *poDataset = (GDALDataset*)GDALOpen(filePath, GA_ReadOnly); 
	if(poDataset==NULL)
	{
		printf("Open file using gdal failed ! \n");
		return 0;
	}

	GDALRasterBand *poBand  = poDataset->GetRasterBand(bandId);
	wd = poBand->GetXSize();
	ht = poBand->GetYSize();
	
	
	
	int ntype = poBand->GetRasterDataType();
	*pBuffer=(unsigned char*)malloc(sizeof(unsigned char)*wd*ht);

	unsigned short* pUShort = NULL;
	short* pShort = NULL;
	unsigned int* pUInt = NULL;
	int* pInt = NULL; 
	float* pFloat = NULL;
	double* pDouble = NULL;
	double dMin, dMax;

	switch(ntype)
	{
	case 1: //byte
		poBand->RasterIO(GF_Read,0,0,wd,ht,*pBuffer,wd,ht,GDT_Byte,0,0);
		break;
	case 2: //16 bit unsigned int
		pUShort = (unsigned short*)malloc(ht*wd*2);
		poBand->RasterIO(GF_Read,0,0,wd,ht,pUShort,wd,ht,GDT_UInt16,0,0);
		GetArrayMinMax(pUShort, wd*ht, dMin, dMax);
		for(int i=0; i<wd*ht; i++)
			(*pBuffer)[i] = (pUShort[i] - dMin) / (dMax-dMin)*255;
		break;
	case 3: //16 bit signed int
		pShort = (short*)malloc(ht*wd*2);
		poBand->RasterIO(GF_Read,0,0,wd,ht,pShort,wd,ht,GDT_Int16,0,0);
		GetArrayMinMax(pShort, wd*ht, dMin, dMax);
		for(int i=0; i<wd*ht; i++)
			(*pBuffer)[i] = (pShort[i] - dMin) / (dMax-dMin)*255;
		break;
	case 4: //32 bit unsigned int
		pUInt = (unsigned int*)malloc(ht*wd*4);
		poBand->RasterIO(GF_Read,0,0,wd,ht,pUInt,wd,ht,GDT_UInt32,0,0);
		GetArrayMinMax(pUInt, wd*ht, dMin, dMax);
		for(int i=0; i<wd*ht; i++)
			(*pBuffer)[i] = (pUInt[i] - dMin) / (dMax-dMin)*255;
		break;
	case 5: //32 bit signed int
		pInt = (int*)malloc(ht*wd*4);
		poBand->RasterIO(GF_Read,0,0,wd,ht,pInt,wd,ht,GDT_Int32,0,0);
		GetArrayMinMax(pInt, wd*ht, dMin, dMax);
		for(int i=0; i<wd*ht; i++)
			(*pBuffer)[i] = (pInt[i] - dMin) / (dMax-dMin)*255;
		break;
	case 6: //32 bit float
		pFloat = (float*)malloc(ht*wd*4);
		poBand->RasterIO(GF_Read,0,0,wd,ht,pUInt,wd,ht,GDT_Float32,0,0);
		GetArrayMinMax(pFloat, wd*ht, dMin, dMax);
		for(int i=0; i<wd*ht; i++)
			(*pBuffer)[i] = (pFloat[i] - dMin) / (dMax-dMin)*255;
		break;
	case 7: //64 bit float
		pDouble = (double*)malloc(ht*wd*8);
		poBand->RasterIO(GF_Read,0,0,wd,ht,pDouble,wd,ht,GDT_Float64,0,0);
		GetArrayMinMax(pDouble, wd*ht, dMin, dMax);
		for(int i=0; i<wd*ht; i++)
			(*pBuffer)[i] = (pDouble[i] - dMin) / (dMax-dMin)*255;
		break;
	default:
		break;		
	}

	free(pUShort);
	free(pShort);
	free(pUInt);
	free(pInt);
	free(pFloat);
	free(pDouble);

	GDALClose( (GDALDatasetH) poDataset );

	return 1;
}

void RegisterAll()
{
	if( GDALGetDriverByName("GTiff") == NULL) 
	{
		GDALAllRegister();
	}
}

/**
所有subset转成一个tif
@param fileName    输入文件名
@param outfile     输出文件名
@param pBandIndex  输出数据集，用数组表示，从1开始，例如[1,3,5]，传入NULL表示全部数据集
@param pBandCount  输出数据集个数，如果pBandIndex设置为NULL，此参数无作用
@return
*/
void subsets2tif(const char* fileName,const char* outfile,int* pBandIndex, int pBandCount)
{
	int i,j;
	char * outFileName=NULL; //输出文件名

	if(outfile == NULL) 
	{
		int charNum = strlen(fileName);
		charNum += 5;
		outFileName = new char[charNum];
		strcpy(outFileName,fileName);
		strcat(outFileName,".tif") ;

	}
	else 
	{
		int charNum = strlen(outfile);
		outFileName = new char[charNum+1];
		strcpy(outFileName,outfile);
	}

	RegisterAll();//注册类型，打开影像必须加入此句
	
	GDALDataset *pDataSet = (GDALDataset *) GDALOpen(fileName, GA_ReadOnly);

	if(pDataSet == NULL) 
	{
		printf("不能打开该文件，请检查文件是否存在！");
		return ;
	}

	//获取子数据集
	char ** papszSUBDATASETS =pDataSet->GetMetadata("SUBDATASETS");
	int subdsNumber =  papszSUBDATASETS == NULL?1:CSLCount(papszSUBDATASETS)/2;
	char** vSubDataSets=new char*[subdsNumber]; //子数据集名称列表
	char** vSubDataDesc=new char*[subdsNumber];  //子数据集描述列表
	char * papszMetadata = NULL;

	//如果没有子数据集，则其本身就是一个数据集
	if(papszSUBDATASETS == NULL) 
	{
		const char* Metadata = GDALGetDriverShortName((GDALDriverH)pDataSet);
		vSubDataSets[0]= new char[strlen(Metadata)+1];
		vSubDataDesc[0]= new char[strlen(Metadata)+1];
		strcpy(vSubDataSets[0], Metadata);
		strcpy(vSubDataDesc[0], Metadata);
	} 
	else 
	{
		int iCount = CSLCount(papszSUBDATASETS);  //计算子数据集的个数
		if(iCount <= 0) 
		{   //没有子数据集,则返回
			GDALClose((GDALDriverH)pDataSet);
			return;
		}

		//将子数据集压入列表
		for(i=0,j=0; papszSUBDATASETS[i] != NULL; i++) 
		{
			
			if(i%2 != 0) 
			{
				continue;
			}
			

			char * setInfo = papszSUBDATASETS[i];
			setInfo = strstr(setInfo,"=");

			if(setInfo[0] == '=') 
			{
				memmove(setInfo, setInfo+1, strlen(setInfo)); //提取元数据中子数据集名称
				vSubDataSets[j]=new char[strlen(setInfo)+1];
				strcpy(vSubDataSets[j],setInfo);
			}

			char * descptn = papszSUBDATASETS[i+1];
			descptn = strstr(descptn,"=");

			if(descptn[0] == '=') 
			{
				memmove(descptn, descptn+1, strlen(descptn)); //提取元数据中子数据集描述
				vSubDataDesc[j]=new char[strlen(descptn)+1];
				strcpy(vSubDataDesc[j],descptn);
			}

			printf("%s \n", vSubDataDesc[j]);

			j++;
		}//end for
	}

	int Width,Height;
	GDALDriver *poDriver;
	const char *pszFormat ="GTiff";
	poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);

	if(poDriver == NULL) 
	{
		return;
	}

	//将每个子数据集转为对应的tif假设所有数据集大小相同
	bool init = false;
	//int iBands = pBandIndex!=NULL?pBandCount:subdsNumber;
	int iBands = 2; //for modis
	int realBandCount=0;

	for(i = 0; i<subdsNumber; i++) 
	{
		if( i!=4 && i!=7 )
			continue;
		
		if(pBandIndex!=NULL) 
		{
			bool runNext =true;

			for(j = 0; j < pBandCount; j++) 
			{
				if(pBandIndex[j]==(i+1)) 
				{
					realBandCount = j+1;
					runNext=false;
					break;
				}
			}

			if(runNext)
			{
				continue;
			}

		}
		else 
		{
			//realBandCount = i+1;
			realBandCount++;
		}

		char* dsName = vSubDataSets[i];
		GDALDataset *subDs = (GDALDataset *) GDALOpen(dsName, GA_ReadOnly);

		if(subDs == NULL) 
		{
			continue;
		}

		Width = subDs->GetRasterXSize();
		Height = subDs->GetRasterYSize();
		int nRasterCount = subDs->GetRasterCount();
		printf("Count: %d \n", nRasterCount);
		void *subData;
		GDALRasterBand *poBand = subDs->GetRasterBand(1);
		GDALDataType eBufType= poBand->GetRasterDataType();

		switch(eBufType) 
		{
		case GDT_Byte:
			subData = new char[Width*Height];
			break;

		case GDT_Int16:
			subData = new short[Width*Height];
			break;

		case GDT_UInt16:
			subData = new unsigned short[Width*Height];
			break;

		case GDT_UInt32:
			subData = new unsigned int[Width*Height];
			break;

		case GDT_Int32:
			subData = new int[Width*Height];
			break;

		case GDT_Float32:
			subData = new float[Width*Height];
			break;

		case GDT_Float64:
			subData = new double[Width*Height];
			break;

		default:
			subData = new double[Width*Height];
			break;
		}

		poBand->RasterIO(GF_Read,0,0,Width,Height,subData,Width,Height,eBufType,0,0);
		
		//写入影像
		GDALDataset *pDstDs = NULL;
		if(!init) 
		{ 
			//创建影像
			char **papszOptions = NULL;
			//设置bsq或者	BIP bsq:BAND,bip:PIXEL
			papszOptions = CSLSetNameValue(papszOptions, "INTERLEAVE", "BAND");  
			pDstDs = poDriver->Create(outFileName,Width,Height,iBands,eBufType,papszOptions);
			double geos[6];
			subDs->GetGeoTransform(geos);//变换参数
			char *pszProjection = NULL;
			char *pszPrettyWkt = NULL;
			//获取ds1坐标
			OGRSpatialReferenceH  hSRS;

			if(GDALGetProjectionRef(subDs) != NULL) 
			{
				pszProjection = (char *) GDALGetProjectionRef(subDs);
				hSRS = OSRNewSpatialReference(NULL);

				if(OSRImportFromWkt(hSRS, &pszProjection) == CE_None) 
				{
					OSRExportToPrettyWkt(hSRS, &pszPrettyWkt, FALSE);
					pDstDs->SetProjection(pszPrettyWkt);
				}
			}

			//设置坐标
			pDstDs->SetGeoTransform(geos);
			CPLFree(pszPrettyWkt);
			pDstDs->SetMetadata(subDs->GetMetadata());
			pDstDs->FlushCache();
			init = true;
		}
		else 
		{
			//打开影像
			pDstDs = (GDALDataset *) GDALOpen(outFileName, GA_Update);
		}

		GDALRasterBand *poBandOut;
		poBandOut = pDstDs->GetRasterBand(realBandCount);
		poBandOut->SetScale(poBand->GetScale());
		poBandOut->SetOffset(poBand->GetOffset());
		poBandOut->SetUnitType(poBand->GetUnitType());
		poBandOut->SetColorInterpretation(poBand->GetColorInterpretation());
		poBandOut->SetDescription(poBand->GetDescription());
		poBandOut->SetNoDataValue(poBand->GetNoDataValue());
		
		//GDT_Float32和 **OutputImg 类型要对应！
		poBandOut->RasterIO(GF_Write,0, 0, Width, Height,subData,Width,Height,eBufType,0,0);
		pDstDs->FlushCache();
		//关闭
		GDALClose(subDs);
		GDALClose(pDstDs);
		delete[] subData;
	}

	GDALClose(pDataSet);

	if(papszMetadata) 
	{
		delete papszMetadata;
	}

	for(i = 0; i < subdsNumber; i++) 
	{
		delete[] vSubDataDesc[i];
		delete[] vSubDataSets[i];
	}

	delete[] vSubDataDesc;
	vSubDataDesc = NULL;
	delete[] vSubDataSets;
	vSubDataSets =NULL;
	delete[] outFileName;

	printf("Finished! \n");
}

int  ReadMod03(char* filePath, double** lon, double** lat, int& ht, int& wd)
{



	return 1;
}

int  ReadSubDataInfo(char* filePath, char** subDataDesc)
{
	int nSubData = 0;

	GDALAllRegister();		

	GDALDataset *pDataSet = (GDALDataset *) GDALOpen(filePath, GA_ReadOnly);
	if(pDataSet == NULL) 
	{
		printf("不能打开该文件，请检查文件是否存在！");
		return 0;
	}

	//获取子数据集
	char ** papszSUBDATASETS =pDataSet->GetMetadata("SUBDATASETS");
	int     subdsNumber =  papszSUBDATASETS == NULL?1:CSLCount(papszSUBDATASETS)/2;
	char**  vSubDataSets=new char*[subdsNumber]; //子数据集名称列表
	char**  vSubDataDesc=new char*[subdsNumber];  //子数据集描述列表
	char *  papszMetadata = NULL;

	//如果没有子数据集，则其本身就是一个数据集
	if(papszSUBDATASETS == NULL) 
	{
		const char* Metadata = GDALGetDriverShortName((GDALDriverH)pDataSet);
		vSubDataSets[0]= new char[strlen(Metadata)+1];
		vSubDataDesc[0]= new char[strlen(Metadata)+1];
		strcpy(vSubDataSets[0], Metadata);
		strcpy(vSubDataDesc[0], Metadata);
	} 
	else 
	{
		int iCount = CSLCount(papszSUBDATASETS);  //计算子数据集的个数
		if(iCount <= 0) 
		{   //没有子数据集,则返回
			GDALClose((GDALDriverH)pDataSet);
			return 0;
		}

		//将子数据集压入列表
		for(int i=0,j=0; papszSUBDATASETS[i] != NULL; i++) 
		{

			if(i%2 != 0) 
			{
				continue;
			}

			char * setInfo = papszSUBDATASETS[i];
			setInfo = strstr(setInfo,"=");
			if(setInfo[0] == '=') 
			{
				memmove(setInfo, setInfo+1, strlen(setInfo)); //提取元数据中子数据集名称
				vSubDataSets[j]=new char[strlen(setInfo)+1];
				strcpy(vSubDataSets[j],setInfo);
			}

			char * descptn = papszSUBDATASETS[i+1];
			descptn = strstr(descptn,"=");
			if(descptn[0] == '=') 
			{
				memmove(descptn, descptn+1, strlen(descptn)); //提取元数据中子数据集描述
				vSubDataDesc[j]=new char[strlen(descptn)+1];
				strcpy(vSubDataDesc[j],descptn);
			}
			printf("%s \n", vSubDataDesc[j]);

			//strcpy( subDataDesc[nSubData], vSubDataDesc[j] );
			strcpy( subDataDesc[nSubData], vSubDataSets[j] );
			nSubData++;

			j++;
		}//end for
	}

	GDALClose(pDataSet);

	if(papszMetadata) 
	{
		delete papszMetadata;
	}
	for(int i=0; i<subdsNumber; i++) 
	{
		delete[] vSubDataDesc[i];
		delete[] vSubDataSets[i];
	}
	delete[] vSubDataDesc;
	vSubDataDesc = NULL;
	delete[] vSubDataSets;
	vSubDataSets =NULL;

	return nSubData;
}

int  ReadSubDataInfo(char* filePath, vector<string>& subDataSets, vector<string>& subDataDesc)
{
	GDALAllRegister();		

	GDALDataset *pDataSet = (GDALDataset *) GDALOpen(filePath, GA_ReadOnly);
	if(pDataSet == NULL) 
	{
		printf("不能打开该文件，请检查文件是否存在！");
		return 0;
	}

	//获取子数据集
	char ** papszSUBDATASETS =pDataSet->GetMetadata("SUBDATASETS");
	//char ** papszSUBDATASETS =pDataSet->GetMetadata("GEOLOCATION");    
	int subdsNumber =  papszSUBDATASETS == NULL?1:CSLCount(papszSUBDATASETS)/2;
	char** vSubDataSets=new char*[subdsNumber]; //子数据集名称列表
	char** vSubDataDesc=new char*[subdsNumber];  //子数据集描述列表
	char * papszMetadata = NULL;

	//如果没有子数据集，则其本身就是一个数据集
	if(papszSUBDATASETS == NULL) 
	{
		const char* Metadata = GDALGetDriverShortName((GDALDriverH)pDataSet);
		vSubDataSets[0]= new char[strlen(Metadata)+1];
		vSubDataDesc[0]= new char[strlen(Metadata)+1];
		strcpy(vSubDataSets[0], Metadata);
		strcpy(vSubDataDesc[0], Metadata);
	} 
	else 
	{
		int iCount = CSLCount(papszSUBDATASETS);  //计算子数据集的个数
		if(iCount <= 0) 
		{   //没有子数据集,则返回
			GDALClose((GDALDriverH)pDataSet);
			return 0;
		}

		//将子数据集压入列表
		for(int i=0,j=0; papszSUBDATASETS[i] != NULL; i++) 
		{

			if(i%2 != 0) 
			{
				continue;
			}

			char * setInfo = papszSUBDATASETS[i];
			setInfo = strstr(setInfo,"=");
			if(setInfo[0] == '=') 
			{
				memmove(setInfo, setInfo+1, strlen(setInfo)); //提取元数据中子数据集名称
				vSubDataSets[j]=new char[strlen(setInfo)+1];
				strcpy(vSubDataSets[j],setInfo);
			}

			char * descptn = papszSUBDATASETS[i+1];
			descptn = strstr(descptn,"=");
			if(descptn[0] == '=') 
			{
				memmove(descptn, descptn+1, strlen(descptn)); //提取元数据中子数据集描述
				vSubDataDesc[j]=new char[strlen(descptn)+1];
				strcpy(vSubDataDesc[j],descptn);
			}
			printf("%s \n", vSubDataDesc[j]);

			
			//string  sfile(vSubDataSets[j]);
			//convert char* to string			
			//subDataSets.push_back("Hello!");
			//string  sfile = vSubDataSets[j];
			//subDataSets.push_back(sfile);
			//sfile = vSubDataDesc[j];
			//subDataDesc.push_back(sfile);			
			//subDataSets.push_back( vSubDataSets[j] );
			//subDataDesc.push_back( vSubDataDesc[j] );
			
			j++;
		}//end for
	}

	GDALClose(pDataSet);

	if(papszMetadata) 
	{
		delete papszMetadata;
	}
	for(int i=0; i<subdsNumber; i++) 
	{
		delete[] vSubDataDesc[i];
		delete[] vSubDataSets[i];
	}
	delete[] vSubDataDesc;
	vSubDataDesc = NULL;
	delete[] vSubDataSets;
	vSubDataSets =NULL;
	
	return 1;
}


int  ReadSubDataBandSum(char* subDataDesc)
{
	GDALAllRegister();

	GDALDataset *subDs = (GDALDataset *) GDALOpen(subDataDesc, GA_ReadOnly);
	if(subDs == NULL) 
	{
		return 0;
	}
	//int Width = subDs->GetRasterXSize();
	//int Height = subDs->GetRasterYSize();
	int nRasterCount = subDs->GetRasterCount();
	//printf("Count: %d \n", nRasterCount);

	GDALClose(subDs);

	return nRasterCount;
}

int  ReadSubDataCaliPara(char* subDataDesc, vector<double>& scale, vector<double>& offset)
{
	GDALAllRegister();

	GDALDataset *subDs = (GDALDataset*)GDALOpen(subDataDesc, GA_ReadOnly);
	if(subDs == NULL) 
	{
		return 0;
	}

	int bandcnt = subDs->GetRasterCount();
	double reflectance_scale;
	double reflectance_offset;

	//提取相应的参数
	char **dataset_str= subDs->GetMetadata();
	for(int j = 0; dataset_str[j] != NULL; j++ )
	{
		
		//printf("%s \n", dataset_str[j]);
		int nLen = strlen(dataset_str[j]);

		char* ts = new char[nLen];
		//const int nLen = 2048;
		//char ts[nLen];
		memset(ts, '\0', nLen);
		strcpy(ts, dataset_str[j]);

		char* pdes = strstr(ts, "=");
		int nCount = pdes-ts;
		
		char* tmpstr = new char[nLen];
		//char tmpstr[nLen];
		memset(tmpstr, '\0', nLen);
		strncpy(tmpstr, ts, nCount);
		
		char* tmp_infostr = new char[nLen];
		//char tmp_infostr[nLen]; 
		memset(tmp_infostr, '\0', nLen);
        strcpy( tmp_infostr, ts+nCount+1 );
		

		//std::string tmp_longstr = dataset_str[j];
		//std::string tmp_infostr = tmp_longstr.substr(tmp_longstr.find_first_of("=")+1);
		//std::string tmpstr = tmp_longstr.substr(0,tmp_longstr.find_first_of("="));

		//if( tmpstr.compare("reflectance_scales")==0)
		if( strcmp(tmpstr,"reflectance_scales")==0 )
		{
			//printf("%s\n",tmpstr.c_str());
			for (int k=0; k<bandcnt; k++)
			{
				char tmp_shortinfostr[256];
				memset(tmp_shortinfostr, '\0', 256);				
				char* pdes = strstr(tmp_infostr, ",");				
				if(pdes!=NULL)
				{
					strncpy(tmp_shortinfostr, tmp_infostr, pdes-tmp_infostr);
					strcpy(tmp_infostr, pdes+1);
				}
				else
				{
					strcpy(tmp_shortinfostr, tmp_infostr);
				}
				sscanf(tmp_shortinfostr,"%lf", &reflectance_scale);
				scale.push_back(reflectance_scale);

				//std::string tmp_shortinfostr = tmp_infostr.substr(0,tmp_infostr.find_first_of(", "));
				//tmp_infostr = tmp_infostr.substr(tmp_infostr.find_first_of(", ")+2);
				//const char *tmpstring=tmp_shortinfostr.c_str();
				//sscanf(tmpstring,"%f", &reflectance_scale);
			}
		}
		
		//if (tmpstr.compare("reflectance_offsets")==0)
		if( strcmp(tmpstr,"reflectance_offsets")==0 )
		{
			//printf("%s\n",tmpstr.c_str());
			for (int k=0; k<bandcnt; k++)
			{
				char tmp_shortinfostr[256];
				memset(tmp_shortinfostr, '\0', 256);
				char* pdes = strstr(tmp_infostr, ",");

				if(pdes!=NULL)
				{
					strncpy(tmp_shortinfostr, tmp_infostr, pdes-tmp_infostr);
					strcpy(tmp_infostr, pdes+1);
				}
				else
				{
					strcpy(tmp_shortinfostr, tmp_infostr);
				}
				sscanf(tmp_shortinfostr,"%lf", &reflectance_offset);
				offset.push_back(reflectance_offset);

				//std::string tmp_shortinfostr = tmp_infostr.substr(0,tmp_infostr.find_first_of(", "));
				//tmp_infostr = tmp_infostr.substr(tmp_infostr.find_first_of(", ")+2);
				//const char *tmpstring=tmp_shortinfostr.c_str();
				//sscanf(tmpstring,"%f", &reflectance_scale);
			}
		}

		//delete ts;
		//delete tmpstr;
		//delete tmp_infostr;
	}

	return 1;
}


int  ReadSubDataFloat(char* subDataDesc, float** pBuffer, int* ht, int* wd, int band)
{
	int i,j;
	//printf("Read subdata: %s  band:%d \n", subDataDesc, band);
	GDALAllRegister();
	GDALDataset *subDs = (GDALDataset *) GDALOpen(subDataDesc, GA_ReadOnly);
	if(subDs == NULL) 
	{
		return 0;
	}

	int Width = subDs->GetRasterXSize();
	int Height = subDs->GetRasterYSize();
	int nRasterCount = subDs->GetRasterCount();
	printf("Count: %d \n", nRasterCount);

	GDALRasterBand *poBand = subDs->GetRasterBand(band);
	GDALDataType eBufType= poBand->GetRasterDataType();

	assert(eBufType==GDT_Float32);

	*pBuffer = ( float* )malloc( Width*Height*sizeof(float) );
	poBand->RasterIO(GF_Read,0,0,Width,Height,*pBuffer,Width,Height,eBufType,0,0);

	*ht = Height;
	*wd = Width;

	GDALClose(subDs);

	return 1;
}

int  ReadSubDataShort(char* subDataDesc, short** pBuffer, int* ht, int* wd, int band)
{
	int i,j;
	//printf("Read subdata: %s  band:%d \n", subDataDesc, band);
	GDALAllRegister();
	GDALDataset *subDs = (GDALDataset *) GDALOpen(subDataDesc, GA_ReadOnly);
	if(subDs == NULL) 
	{
		return 0;
	}

	int Width = subDs->GetRasterXSize();
	int Height = subDs->GetRasterYSize();
	int nRasterCount = subDs->GetRasterCount();
	printf("Count: %d \n", nRasterCount);

	GDALRasterBand *poBand = subDs->GetRasterBand(band);
	GDALDataType eBufType= poBand->GetRasterDataType();

	assert(eBufType==GDT_Int16);

	*pBuffer = ( short* )malloc( Width*Height*sizeof(short) );
	poBand->RasterIO(GF_Read,0,0,Width,Height,*pBuffer,Width,Height,eBufType,0,0);

	*ht = Height;
	*wd = Width;

	GDALClose(subDs);

	return 1;
}


int  ReadSubDataUShort(char* subDataDesc, unsigned short** pBuffer, int* ht, int* wd, int band)
{
	int i,j;
	//printf("Read subdata: %s  band:%d \n", subDataDesc, band);
	GDALAllRegister();
	GDALDataset *subDs = (GDALDataset *) GDALOpen(subDataDesc, GA_ReadOnly);
	if(subDs == NULL) 
	{
		return 0;
	}

	int Width = subDs->GetRasterXSize();
	int Height = subDs->GetRasterYSize();
	int nRasterCount = subDs->GetRasterCount();
	printf("Count: %d \n", nRasterCount);
	
	GDALRasterBand *poBand = subDs->GetRasterBand(band);
	GDALDataType eBufType= poBand->GetRasterDataType();

	assert(eBufType==GDT_UInt16);

	*pBuffer = ( unsigned short* )malloc( Width*Height*sizeof(unsigned short) );
	poBand->RasterIO(GF_Read,0,0,Width,Height,*pBuffer,Width,Height,eBufType,0,0);
    
	*ht = Height;
	*wd = Width;

	GDALClose(subDs);

	return 1;
}


int  ReadSubData(char* subDataDesc, double** pBuffer, int* ht, int* wd,int band)
{
	int i,j;

	printf("Read subdata: %s  band:%d \n", subDataDesc, band);

	GDALAllRegister();

	GDALDataset *subDs = (GDALDataset *) GDALOpen(subDataDesc, GA_ReadOnly);
	if(subDs == NULL) 
	{
		return 0;
	}

	int Width = subDs->GetRasterXSize();
	int Height = subDs->GetRasterYSize();
	int nRasterCount = subDs->GetRasterCount();
	printf("Count: %d \n", nRasterCount);
	void *subData;
	GDALRasterBand *poBand = subDs->GetRasterBand(band);
	GDALDataType eBufType= poBand->GetRasterDataType();

	switch(eBufType) 
	{
	case GDT_Byte:
		subData = new char[Width*Height];
		break;

	case GDT_Int16:
		subData = new short[Width*Height];
		break;

	case GDT_UInt16:
		subData = new unsigned short[Width*Height];
		break;

	case GDT_UInt32:
		subData = new unsigned int[Width*Height];
		break;

	case GDT_Int32:
		subData = new int[Width*Height];
		break;

	case GDT_Float32:
		subData = new float[Width*Height];
		break;

	case GDT_Float64:
		subData = new double[Width*Height];
		break;

	default:
		subData = new double[Width*Height];
		break;
	}

	poBand->RasterIO(GF_Read,0,0,Width,Height,subData,Width,Height,eBufType,0,0);
	
	*ht = Height;
	*wd = Width;
	*pBuffer = (double*)malloc(Width*Height*sizeof(double));

	//int floatSize = sizeof(float);

	switch(eBufType) 
	{
	case GDT_Byte:
		for( j=0; j<Height; j++)
			for( i=0; i<Width; i++)
			{
				(*pBuffer)[j*Width+i] = *( (unsigned char*)(subData)+(j*Width+i) ); 
			}
		break;

	case GDT_Int16:
		for( j=0; j<Height; j++)
			for( i=0; i<Width; i++)
			{
				(*pBuffer)[j*Width+i] = *( (short*)(subData)+(j*Width+i) ); 
			}
		break;

	case GDT_UInt16:
		for( j=0; j<Height; j++)
			for( i=0; i<Width; i++)
			{
				(*pBuffer)[j*Width+i] = *( (unsigned short*)(subData)+(j*Width+i) ); 
			}
		break;

	case GDT_UInt32:
		for( j=0; j<Height; j++)
			for( i=0; i<Width; i++)
			{
				(*pBuffer)[j*Width+i] = *( (unsigned int*)(subData)+(j*Width+i) ); 
			}
		break;

	case GDT_Int32:
		for( j=0; j<Height; j++)
			for( i=0; i<Width; i++)
			{
				(*pBuffer)[j*Width+i] = *( (int*)(subData)+(j*Width+i)); 
			}
		break;

	case GDT_Float32:
		for( j=0; j<Height; j++)
			for( i=0; i<Width; i++)
			{
				(*pBuffer)[j*Width+i] = *( (float*)(subData)+(j*Width+i) ); 
			}
		break;

	case GDT_Float64:
		for( j=0; j<Height; j++)
			for( i=0; i<Width; i++)
			{
				(*pBuffer)[j*Width+i] = *( (double*)(subData)+(j*Width+i) ); 
			}
		break;

	default:
		for( j=0; j<Height; j++)
			for( i=0; i<Width; i++)
			{
				(*pBuffer)[j*Width+i] = *( (double*)(subData)+(j*Width+i)*sizeof(double) ); 
			}
		break;
	}


	GDALClose(subDs);
	delete[] subData;

	return 1;
}



/************************************************************************/
/*                        GDALInfoReportCorner()                        */
/************************************************************************/

static int 
GDALInfoReportCorner( GDALDatasetH hDataset, 
					 OGRCoordinateTransformationH hTransform,
					 const char * corner_name,
					 double x, double y )

{
	double	dfGeoX, dfGeoY;
	double	adfGeoTransform[6];

	printf( "%-11s ", corner_name );

	/* -------------------------------------------------------------------- */
	/*      Transform the point into georeferenced coordinates.             */
	/* -------------------------------------------------------------------- */
	if( GDALGetGeoTransform( hDataset, adfGeoTransform ) == CE_None )
	{
		dfGeoX = adfGeoTransform[0] + adfGeoTransform[1] * x
			+ adfGeoTransform[2] * y;
		dfGeoY = adfGeoTransform[3] + adfGeoTransform[4] * x
			+ adfGeoTransform[5] * y;
	}

	else
	{
		printf( "(%7.1f,%7.1f)\n", x, y );
		return FALSE;
	}

	/* -------------------------------------------------------------------- */
	/*      Report the georeferenced coordinates.                           */
	/* -------------------------------------------------------------------- */
	if( ABS(dfGeoX) < 181 && ABS(dfGeoY) < 91 )
	{
		printf( "(%12.7f,%12.7f) ", dfGeoX, dfGeoY );

	}
	else
	{
		printf( "(%12.3f,%12.3f) ", dfGeoX, dfGeoY );
	}

	/* -------------------------------------------------------------------- */
	/*      Transform to latlong and report.                                */
	/* -------------------------------------------------------------------- */
	if( hTransform != NULL 
		&& OCTTransform(hTransform,1,&dfGeoX,&dfGeoY,NULL) )
	{

		printf( "(%s,", GDALDecToDMS( dfGeoX, "Long", 2 ) );
		printf( "%s)", GDALDecToDMS( dfGeoY, "Lat", 2 ) );
	}

	printf( "\n" );

	return TRUE;
}

int  GdalInfo(char* pszFilename)
{
	GDALDatasetH		hDataset;
	GDALRasterBandH		hBand;
	int					i, iBand;
	double				adfGeoTransform[6];
	GDALDriverH			hDriver;
	char				**papszMetadata;
	int                 bComputeMinMax = FALSE, bSample = FALSE;
	int                 bShowGCPs = TRUE, bShowMetadata = TRUE, bShowRAT=TRUE;
	int                 bStats = FALSE, bApproxStats = TRUE, iMDD;
	int                 bShowColorTable = TRUE, bComputeChecksum = FALSE;
	int                 bReportHistograms = FALSE;
	int                 bReportProj4 = FALSE;
	int                 nSubdataset = -1;
	char				**papszExtraMDDomains = NULL, **papszFileList;
	const char			*pszProjection = NULL;
	OGRCoordinateTransformationH hTransform = NULL;
	int					bShowFileList = TRUE;

	/* Check that we are running against at least GDAL 1.5 */
	/* Note to developers : if we use newer API, please change the requirement */
	
	if (atoi(GDALVersionInfo("VERSION_NUM")) < 1500)
	{
		printf("At least, GDAL >= 1.5.0 is required for this version, "
			"which was compiled against GDAL %s\n",  GDAL_RELEASE_NAME);
		
		return 0; 
		//exit(1);
	}
	
	GDALAllRegister();

	/* -------------------------------------------------------------------- */
	/*      Open dataset.                                                   */
	/* -------------------------------------------------------------------- */
	hDataset = GDALOpen( pszFilename, GA_ReadOnly );

	if( hDataset == NULL )
	{
		fprintf( stderr,
			"gdalinfo failed - unable to open '%s'.\n",
			pszFilename );

		//CSLDestroy( argv );
		//CSLDestroy( papszExtraMDDomains );

		GDALDumpOpenDatasets( stderr );

		GDALDestroyDriverManager();

		CPLDumpSharedList( NULL );

		exit( 1 );
	}

	/* -------------------------------------------------------------------- */
	/*      Read specified subdataset if requested.                         */
	/* -------------------------------------------------------------------- */
	if ( nSubdataset > 0 )
	{
		char **papszSubdatasets = GDALGetMetadata( hDataset, "SUBDATASETS" );
		int nSubdatasets = CSLCount( papszSubdatasets );

		if ( nSubdatasets > 0 && nSubdataset <= nSubdatasets )
		{
			char szKeyName[1024];
			char *pszSubdatasetName;

			snprintf( szKeyName, sizeof(szKeyName),
				"SUBDATASET_%d_NAME", nSubdataset );
			szKeyName[sizeof(szKeyName) - 1] = '\0';
			pszSubdatasetName =
				CPLStrdup( CSLFetchNameValue( papszSubdatasets, szKeyName ) );
			GDALClose( hDataset );
			hDataset = GDALOpen( pszSubdatasetName, GA_ReadOnly );
			CPLFree( pszSubdatasetName );
		}
		else
		{
			fprintf( stderr,
				"gdalinfo warning: subdataset %d of %d requested. "
				"Reading the main dataset.\n",
				nSubdataset, nSubdatasets );

		}
	}

	/* -------------------------------------------------------------------- */
	/*      Report general info.                                            */
	/* -------------------------------------------------------------------- */
	hDriver = GDALGetDatasetDriver( hDataset );
	printf( "Driver: %s/%s\n",
		GDALGetDriverShortName( hDriver ),
		GDALGetDriverLongName( hDriver ) );

	/*
	papszFileList = GDALGetFileList( hDataset );
	if( CSLCount(papszFileList) == 0 )
	{
		printf( "Files: none associated\n" );
	}
	else
	{
		printf( "Files: %s\n", papszFileList[0] );
		if( bShowFileList )
		{
			for( i = 1; papszFileList[i] != NULL; i++ )
				printf( "       %s\n", papszFileList[i] );
		}
	}
	CSLDestroy( papszFileList );
	*/

	printf( "Size is %d, %d\n",
		GDALGetRasterXSize( hDataset ), 
		GDALGetRasterYSize( hDataset ) );

	/* -------------------------------------------------------------------- */
	/*      Report projection.                                              */
	/* -------------------------------------------------------------------- */
	if( GDALGetProjectionRef( hDataset ) != NULL )
	{
		OGRSpatialReferenceH  hSRS;
		char		      *pszProjection;

		pszProjection = (char *) GDALGetProjectionRef( hDataset );

		hSRS = OSRNewSpatialReference(NULL);
		if( OSRImportFromWkt( hSRS, &pszProjection ) == CE_None )
		{
			char	*pszPrettyWkt = NULL;

			OSRExportToPrettyWkt( hSRS, &pszPrettyWkt, FALSE );
			printf( "Coordinate System is:\n%s\n", pszPrettyWkt );
			CPLFree( pszPrettyWkt );
		}
		else
			printf( "Coordinate System is `%s'\n",
			GDALGetProjectionRef( hDataset ) );

		if ( bReportProj4 ) 
		{
			char *pszProj4 = NULL;
			OSRExportToProj4( hSRS, &pszProj4 );
			printf("PROJ.4 string is:\n\'%s\'\n",pszProj4);
			CPLFree( pszProj4 ); 
		}

		OSRDestroySpatialReference( hSRS );
	}

	/* -------------------------------------------------------------------- */
	/*      Report Geotransform.                                            */
	/* -------------------------------------------------------------------- */
	if( GDALGetGeoTransform( hDataset, adfGeoTransform ) == CE_None )
	{
		if( adfGeoTransform[2] == 0.0 && adfGeoTransform[4] == 0.0 )
		{
			printf( "Origin = (%.15f,%.15f)\n",
				adfGeoTransform[0], adfGeoTransform[3] );

			printf( "Pixel Size = (%.15f,%.15f)\n",
				adfGeoTransform[1], adfGeoTransform[5] );
		}
		else
			printf( "GeoTransform =\n"
			"  %.16g, %.16g, %.16g\n"
			"  %.16g, %.16g, %.16g\n", 
			adfGeoTransform[0],
			adfGeoTransform[1],
			adfGeoTransform[2],
			adfGeoTransform[3],
			adfGeoTransform[4],
			adfGeoTransform[5] );
	}

	/* -------------------------------------------------------------------- */
	/*      Report GCPs.                                                    */
	/* -------------------------------------------------------------------- */
	if( bShowGCPs && GDALGetGCPCount( hDataset ) > 0 )
	{
		if (GDALGetGCPProjection(hDataset) != NULL)
		{
			OGRSpatialReferenceH  hSRS;
			char		      *pszProjection;

			pszProjection = (char *) GDALGetGCPProjection( hDataset );

			hSRS = OSRNewSpatialReference(NULL);
			if( OSRImportFromWkt( hSRS, &pszProjection ) == CE_None )
			{
				char	*pszPrettyWkt = NULL;

				OSRExportToPrettyWkt( hSRS, &pszPrettyWkt, FALSE );
				printf( "GCP Projection = \n%s\n", pszPrettyWkt );
				CPLFree( pszPrettyWkt );
			}
			else
				printf( "GCP Projection = %s\n",
				GDALGetGCPProjection( hDataset ) );

			OSRDestroySpatialReference( hSRS );
		}

		for( i = 0; i < GDALGetGCPCount(hDataset); i++ )
		{
			const GDAL_GCP	*psGCP;

			psGCP = GDALGetGCPs( hDataset ) + i;

			printf( "GCP[%3d]: Id=%s, Info=%s\n"
				"          (%.15g,%.15g) -> (%.15g,%.15g,%.15g)\n", 
				i, psGCP->pszId, psGCP->pszInfo, 
				psGCP->dfGCPPixel, psGCP->dfGCPLine, 
				psGCP->dfGCPX, psGCP->dfGCPY, psGCP->dfGCPZ );
		}
	}

	/* -------------------------------------------------------------------- */
	/*      Report metadata.                                                */
	/* -------------------------------------------------------------------- */
	papszMetadata = (bShowMetadata) ? GDALGetMetadata( hDataset, NULL ) : NULL;
	if( bShowMetadata && CSLCount(papszMetadata) > 0 )
	{
		printf( "Metadata:\n" );
		for( i = 0; papszMetadata[i] != NULL; i++ )
		{
			printf( "  %s\n", papszMetadata[i] );
		}
	}

	for( iMDD = 0; bShowMetadata && iMDD < CSLCount(papszExtraMDDomains); iMDD++ )
	{
		papszMetadata = GDALGetMetadata( hDataset, papszExtraMDDomains[iMDD] );
		if( CSLCount(papszMetadata) > 0 )
		{
			printf( "Metadata (%s):\n", papszExtraMDDomains[iMDD]);
			for( i = 0; papszMetadata[i] != NULL; i++ )
			{
				if (EQUALN(papszExtraMDDomains[iMDD], "xml:", 4))
					printf( "%s\n", papszMetadata[i] );
				else
					printf( "  %s\n", papszMetadata[i] );
			}
		}
	}

	/* -------------------------------------------------------------------- */
	/*      Report "IMAGE_STRUCTURE" metadata.                              */
	/* -------------------------------------------------------------------- */
	papszMetadata = (bShowMetadata) ? GDALGetMetadata( hDataset, "IMAGE_STRUCTURE" ) : NULL;
	if( bShowMetadata && CSLCount(papszMetadata) > 0 )
	{
		printf( "Image Structure Metadata:\n" );
		for( i = 0; papszMetadata[i] != NULL; i++ )
		{
			printf( "  %s\n", papszMetadata[i] );
		}
	}

	/* -------------------------------------------------------------------- */
	/*      Report subdatasets.                                             */
	/* -------------------------------------------------------------------- */
	papszMetadata = GDALGetMetadata( hDataset, "SUBDATASETS" );
	if( CSLCount(papszMetadata) > 0 )
	{
		printf( "Subdatasets:\n" );
		for( i = 0; papszMetadata[i] != NULL; i++ )
		{
			printf( "  %s\n", papszMetadata[i] );
		}
	}

	/* -------------------------------------------------------------------- */
	/*      Report geolocation.                                             */
	/* -------------------------------------------------------------------- */
	papszMetadata = (bShowMetadata) ? GDALGetMetadata( hDataset, "GEOLOCATION" ) : NULL;
	if( bShowMetadata && CSLCount(papszMetadata) > 0 )
	{
		printf( "Geolocation:\n" );
		for( i = 0; papszMetadata[i] != NULL; i++ )
		{
			printf( "  %s\n", papszMetadata[i] );
		}
	}

	/* -------------------------------------------------------------------- */
	/*      Report RPCs                                                     */
	/* -------------------------------------------------------------------- */
	papszMetadata = (bShowMetadata) ? GDALGetMetadata( hDataset, "RPC" ) : NULL;
	if( bShowMetadata && CSLCount(papszMetadata) > 0 )
	{
		printf( "RPC Metadata:\n" );
		for( i = 0; papszMetadata[i] != NULL; i++ )
		{
			printf( "  %s\n", papszMetadata[i] );
		}
	}

	/* -------------------------------------------------------------------- */
	/*      Setup projected to lat/long transform if appropriate.           */
	/* -------------------------------------------------------------------- */
	if( GDALGetGeoTransform( hDataset, adfGeoTransform ) == CE_None )
		pszProjection = GDALGetProjectionRef(hDataset);

	if( pszProjection != NULL && strlen(pszProjection) > 0 )
	{
		OGRSpatialReferenceH hProj, hLatLong = NULL;

		hProj = OSRNewSpatialReference( pszProjection );
		if( hProj != NULL )
			hLatLong = OSRCloneGeogCS( hProj );

		if( hLatLong != NULL )
		{
			CPLPushErrorHandler( CPLQuietErrorHandler );
			hTransform = OCTNewCoordinateTransformation( hProj, hLatLong );
			CPLPopErrorHandler();

			OSRDestroySpatialReference( hLatLong );
		}

		if( hProj != NULL )
			OSRDestroySpatialReference( hProj );
	}

	/* -------------------------------------------------------------------- */
	/*      Report corners.                                                 */
	/* -------------------------------------------------------------------- */
	printf( "Corner Coordinates:\n" );
	GDALInfoReportCorner( hDataset, hTransform, "Upper Left", 
		0.0, 0.0 );
	GDALInfoReportCorner( hDataset, hTransform, "Lower Left", 
		0.0, GDALGetRasterYSize(hDataset));
	GDALInfoReportCorner( hDataset, hTransform, "Upper Right", 
		GDALGetRasterXSize(hDataset), 0.0 );
	GDALInfoReportCorner( hDataset, hTransform, "Lower Right", 
		GDALGetRasterXSize(hDataset), 
		GDALGetRasterYSize(hDataset) );
	GDALInfoReportCorner( hDataset, hTransform, "Center", 
		GDALGetRasterXSize(hDataset)/2.0, 
		GDALGetRasterYSize(hDataset)/2.0 );

	if( hTransform != NULL )
	{
		OCTDestroyCoordinateTransformation( hTransform );
		hTransform = NULL;
	}

	/* ==================================================================== */
	/*      Loop over bands.                                                */
	/* ==================================================================== */
	for( iBand = 0; iBand < GDALGetRasterCount( hDataset ); iBand++ )
	{
		double      dfMin, dfMax, adfCMinMax[2], dfNoData;
		int         bGotMin, bGotMax, bGotNodata, bSuccess;
		int         nBlockXSize, nBlockYSize, nMaskFlags;
		double      dfMean, dfStdDev;
		GDALColorTableH	hTable;
		CPLErr      eErr;

		hBand = GDALGetRasterBand( hDataset, iBand+1 );

		if( bSample )
		{
			float afSample[10000];
			int   nCount;

			nCount = GDALGetRandomRasterSample( hBand, 10000, afSample );
			printf( "Got %d samples.\n", nCount );
		}

		GDALGetBlockSize( hBand, &nBlockXSize, &nBlockYSize );
		printf( "Band %d Block=%dx%d Type=%s, ColorInterp=%s\n", iBand+1,
			nBlockXSize, nBlockYSize,
			GDALGetDataTypeName(
			GDALGetRasterDataType(hBand)),
			GDALGetColorInterpretationName(
			GDALGetRasterColorInterpretation(hBand)) );

		if( GDALGetDescription( hBand ) != NULL 
			&& strlen(GDALGetDescription( hBand )) > 0 )
			printf( "  Description = %s\n", GDALGetDescription(hBand) );

		dfMin = GDALGetRasterMinimum( hBand, &bGotMin );
		dfMax = GDALGetRasterMaximum( hBand, &bGotMax );
		if( bGotMin || bGotMax || bComputeMinMax )
		{
			printf( "  " );
			if( bGotMin )
				printf( "Min=%.3f ", dfMin );
			if( bGotMax )
				printf( "Max=%.3f ", dfMax );

			if( bComputeMinMax )
			{
				CPLErrorReset();
				GDALComputeRasterMinMax( hBand, FALSE, adfCMinMax );
				if (CPLGetLastErrorType() == CE_None)
				{
					printf( "  Computed Min/Max=%.3f,%.3f", 
						adfCMinMax[0], adfCMinMax[1] );
				}
			}

			printf( "\n" );
		}

		eErr = GDALGetRasterStatistics( hBand, bApproxStats, bStats, 
			&dfMin, &dfMax, &dfMean, &dfStdDev );
		if( eErr == CE_None )
		{
			printf( "  Minimum=%.3f, Maximum=%.3f, Mean=%.3f, StdDev=%.3f\n",
				dfMin, dfMax, dfMean, dfStdDev );
		}

		if( bReportHistograms )
		{
			int nBucketCount, *panHistogram = NULL;

			eErr = GDALGetDefaultHistogram( hBand, &dfMin, &dfMax, 
				&nBucketCount, &panHistogram, 
				TRUE, GDALTermProgress, NULL );
			if( eErr == CE_None )
			{
				int iBucket;

				printf( "  %d buckets from %g to %g:\n  ",
					nBucketCount, dfMin, dfMax );
				for( iBucket = 0; iBucket < nBucketCount; iBucket++ )
					printf( "%d ", panHistogram[iBucket] );
				printf( "\n" );
				CPLFree( panHistogram );
			}
		}

		if ( bComputeChecksum)
		{
			printf( "  Checksum=%d\n",
				GDALChecksumImage(hBand, 0, 0,
				GDALGetRasterXSize(hDataset),
				GDALGetRasterYSize(hDataset)));
		}

		dfNoData = GDALGetRasterNoDataValue( hBand, &bGotNodata );
		if( bGotNodata )
		{
			if (CPLIsNan(dfNoData))
				printf( "  NoData Value=nan\n" );
			else
				printf( "  NoData Value=%.18g\n", dfNoData );
		}

		if( GDALGetOverviewCount(hBand) > 0 )
		{
			int		iOverview;

			printf( "  Overviews: " );
			for( iOverview = 0; 
				iOverview < GDALGetOverviewCount(hBand);
				iOverview++ )
			{
				GDALRasterBandH	hOverview;
				const char *pszResampling = NULL;

				if( iOverview != 0 )
					printf( ", " );

				hOverview = GDALGetOverview( hBand, iOverview );
				if (hOverview != NULL)
				{
					printf( "%dx%d", 
						GDALGetRasterBandXSize( hOverview ),
						GDALGetRasterBandYSize( hOverview ) );

					pszResampling = 
						GDALGetMetadataItem( hOverview, "RESAMPLING", "" );

					if( pszResampling != NULL 
						&& EQUALN(pszResampling,"AVERAGE_BIT2",12) )
						printf( "*" );
				}
				else
					printf( "(null)" );
			}
			printf( "\n" );

			if ( bComputeChecksum)
			{
				printf( "  Overviews checksum: " );
				for( iOverview = 0; 
					iOverview < GDALGetOverviewCount(hBand);
					iOverview++ )
				{
					GDALRasterBandH	hOverview;

					if( iOverview != 0 )
						printf( ", " );

					hOverview = GDALGetOverview( hBand, iOverview );
					if (hOverview)
						printf( "%d",
						GDALChecksumImage(hOverview, 0, 0,
						GDALGetRasterBandXSize(hOverview),
						GDALGetRasterBandYSize(hOverview)));
					else
						printf( "(null)" );
				}
				printf( "\n" );
			}
		}

		if( GDALHasArbitraryOverviews( hBand ) )
		{
			printf( "  Overviews: arbitrary\n" );
		}

		/*
		nMaskFlags = GDALGetMaskFlags( hBand );
		if( (nMaskFlags & (GMF_NODATA|GMF_ALL_VALID)) == 0 )
		{
			GDALRasterBandH hMaskBand = GDALGetMaskBand(hBand) ;

			printf( "  Mask Flags: " );
			if( nMaskFlags & GMF_PER_DATASET )
				printf( "PER_DATASET " );
			if( nMaskFlags & GMF_ALPHA )
				printf( "ALPHA " );
			if( nMaskFlags & GMF_NODATA )
				printf( "NODATA " );
			if( nMaskFlags & GMF_ALL_VALID )
				printf( "ALL_VALID " );
			printf( "\n" );

			if( hMaskBand != NULL &&
				GDALGetOverviewCount(hMaskBand) > 0 )
			{
				int		iOverview;

				printf( "  Overviews of mask band: " );
				for( iOverview = 0; 
					iOverview < GDALGetOverviewCount(hMaskBand);
					iOverview++ )
				{
					GDALRasterBandH	hOverview;

					if( iOverview != 0 )
						printf( ", " );

					hOverview = GDALGetOverview( hMaskBand, iOverview );
					printf( "%dx%d", 
						GDALGetRasterBandXSize( hOverview ),
						GDALGetRasterBandYSize( hOverview ) );
				}
				printf( "\n" );
			}
		}
		*/

		if( strlen(GDALGetRasterUnitType(hBand)) > 0 )
		{
			printf( "  Unit Type: %s\n", GDALGetRasterUnitType(hBand) );
		}

		if( GDALGetRasterCategoryNames(hBand) != NULL )
		{
			char **papszCategories = GDALGetRasterCategoryNames(hBand);
			int i;

			printf( "  Categories:\n" );
			for( i = 0; papszCategories[i] != NULL; i++ )
				printf( "    %3d: %s\n", i, papszCategories[i] );
		}

		if( GDALGetRasterScale( hBand, &bSuccess ) != 1.0 
			|| GDALGetRasterOffset( hBand, &bSuccess ) != 0.0 )
			printf( "  Offset: %.15g,   Scale:%.15g\n",
			GDALGetRasterOffset( hBand, &bSuccess ),
			GDALGetRasterScale( hBand, &bSuccess ) );

		papszMetadata = (bShowMetadata) ? GDALGetMetadata( hBand, NULL ) : NULL;
		if( bShowMetadata && CSLCount(papszMetadata) > 0 )
		{
			printf( "  Metadata:\n" );
			for( i = 0; papszMetadata[i] != NULL; i++ )
			{
				printf( "    %s\n", papszMetadata[i] );
			}
		}

		papszMetadata = (bShowMetadata) ? GDALGetMetadata( hBand, "IMAGE_STRUCTURE" ) : NULL;
		if( bShowMetadata && CSLCount(papszMetadata) > 0 )
		{
			printf( "  Image Structure Metadata:\n" );
			for( i = 0; papszMetadata[i] != NULL; i++ )
			{
				printf( "    %s\n", papszMetadata[i] );
			}
		}

		if( GDALGetRasterColorInterpretation(hBand) == GCI_PaletteIndex 
			&& (hTable = GDALGetRasterColorTable( hBand )) != NULL )
		{
			int			i;

			printf( "  Color Table (%s with %d entries)\n", 
				GDALGetPaletteInterpretationName(
				GDALGetPaletteInterpretation( hTable )), 
				GDALGetColorEntryCount( hTable ) );

			if (bShowColorTable)
			{
				for( i = 0; i < GDALGetColorEntryCount( hTable ); i++ )
				{
					GDALColorEntry	sEntry;

					GDALGetColorEntryAsRGB( hTable, i, &sEntry );
					printf( "  %3d: %d,%d,%d,%d\n", 
						i, 
						sEntry.c1,
						sEntry.c2,
						sEntry.c3,
						sEntry.c4 );
				}
			}
		}

		if( bShowRAT && GDALGetDefaultRAT( hBand ) != NULL )
		{
			GDALRasterAttributeTableH hRAT = GDALGetDefaultRAT( hBand );

			GDALRATDumpReadable( hRAT, NULL );
		}
	}

	GDALClose( hDataset );

	//CSLDestroy( papszExtraMDDomains );
	//CSLDestroy( argv );

	GDALDumpOpenDatasets( stderr );

	GDALDestroyDriverManager();

	CPLDumpSharedList( NULL );
	CPLCleanupTLS();

	return 1;
}
