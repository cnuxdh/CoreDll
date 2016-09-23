

#ifndef GEO_TIFF_H
#define GEO_TIFF_H

#include "ogr_spatialref.h"


#include <vector>
#include <string>
using namespace std;


#include "exports.h"
#include "gdal_priv.h"
//#include "ogr_spatialref.h"


//#define DLL_EXPORT  _declspec(dllexport)

typedef struct struct_GeoInfo
{
	double left,top;        //ground coordinate of top-left corner 
	double minlon, maxlat;  //(lon,lat) of top-left corner
	double dx,dy;           //resolution;
	int    wd,ht;           //dimension of image
	int    nband;
	int    zoneNumber;      //zone number 
	int    type;            //data type
	const char* projectRef; //
}stGeoInfo;

DLL_EXPORT int ReadGeoFile(char* filename, stGeoInfo& geoData);

//from gdalinfo.c 
DLL_EXPORT int  GdalInfo(char* filePath);

//print information from image using Gdal
DLL_EXPORT int  PrintImageInfo(char* pszFilename);

//retrieve the geo-information from the image file
DLL_EXPORT int  GetGeoInformation(const char* pszFilename, stGeoInfo& geoInfo);

//
DLL_EXPORT GDALDataType  GetDataType(const char* pszFilename);



/*************************   for GEOTiff  ********************************/
//read the band from the file of BYTE type
DLL_EXPORT int  ReadGeoFileByte(char* filePath, int bandId, unsigned char** pBuffer, int& ht, int& wd);
DLL_EXPORT int  ReadGeoFileByte(char* filePath, int bandId, double ratio, unsigned char** pBuffer, int& ht, int& wd);

DLL_EXPORT int  ReadGeoFileShort(char* filePath, int bandId, short** pBuffer, int& ht, int& wd);
//read the band from the file of unsigned short type
DLL_EXPORT int  ReadGeoFileUShort(char* filePath, int bandId, unsigned short** pBuffer, int& ht, int& wd);
//read the band from the file of float type
DLL_EXPORT int  ReadGeoFileFloat(char* filePath, int bandId, float** pBuffer, int& ht, int& wd);
//read the band from the file of int type
DLL_EXPORT int  ReadGeoFileInt(char* filePath, int bandId, int** pBuffer, int& ht, int& wd);
DLL_EXPORT int  ReadGeoFileUInt(char* filePath, int bandId, unsigned int** pBuffer, int& ht, int& wd);

//template funcion
DLL_EXPORT int  ReadGeoFileGeneral(char* filepath, int bandId, void* pBuffer, GDALDataType nType);



//write single band buffer into the file
DLL_EXPORT int  GdalWriteImageByte(char* savePath, unsigned char* pBuffer, int ht, int wd);
DLL_EXPORT int  GdalWriteByte(char* savePath, unsigned char* pBuffer, int ht, int wd, stGeoInfo geoInfo);
DLL_EXPORT int  GdalWriteInt(char* savePath, int* pBuffer, int ht, int wd, stGeoInfo geoInfo);
DLL_EXPORT int  GdalWriteUInt(char* savePath, unsigned int* pBuffer, int ht, int wd, stGeoInfo geoInfo);
DLL_EXPORT int  GdalWriteUShort(char* savePath, unsigned short* pBuffer, int ht, int wd, stGeoInfo geoInfo);



//write r,g,b chanel to the color Tiff file
DLL_EXPORT int  GdalWriteImageByteColor(char* savePath, unsigned char* r, unsigned char* g, unsigned char* b, int ht, int wd);
//write r,g,b chanel to the color Tiff file with Geoinfo
DLL_EXPORT int  GdalWriteImageColor(char* savePath, unsigned char* r, unsigned char* g, unsigned char* b, int ht, int wd,
										stGeoInfo geoInfo);

DLL_EXPORT int  GdalWriteColorImageUShort(char* savePath, 
				unsigned short* r, unsigned short* g, unsigned short* b, int ht, int wd,
				stGeoInfo geoInfo);


//write the float buffer into the GeoTiff file
DLL_EXPORT int  GdalWriteFloat(char* savePath, float* pBuffer, int ht, int wd, stGeoInfo geoInfo);
DLL_EXPORT int  GdalWriteImageUShort(char* savePath, unsigned short* pBuffer, int ht, int wd);
DLL_EXPORT int  GdalWriteJpgCopy(char* savePath, unsigned char* pBuffer, int ht, int wd);
DLL_EXPORT int  GdalWriteFloatLL(char* savePath, float* pBuffer, int ht, int wd, 
								 double minlon, double maxlon, double minlax, double maxlax);

//write r,g,b channels into geotiff file
DLL_EXPORT int  GdalWriteFloatLL(char* savePath, 
								 unsigned char* r, unsigned char* g, unsigned char* b,
								 int ht, int wd, 
								 double minlon, double maxlon, double minlax, double maxlax);


DLL_EXPORT void GroundToLL(double gx, double gy, double& lon, double& lat, stGeoInfo geoInfo);
DLL_EXPORT void GeoImageMosaic(char** files, int nfile);
DLL_EXPORT int  GeotiffMedianFilter(char* filename);
/*************************   for GEOTif  ********************************/



/*************************   for HDF  ********************************/

DLL_EXPORT void subsets2tif(const char* fileName,const char* outfile,int* pBandIndex, int pBandCount);

//read MOD03 product
DLL_EXPORT int  ReadMod03(char* filePath, double** lon, double** lat, int& ht, int& wd);


DLL_EXPORT int  Mod02ToTif(char* mod02File, char* mod03File, char* outDir);


DLL_EXPORT int  ReadSubDataInfo(char* filePath, vector<string>& subDataSets, vector<string>& subDataDesc);

//retrieve the subdata information from the file
DLL_EXPORT int  ReadSubDataInfo(char* filePath, char** subDataDesc);

DLL_EXPORT int  ReadSubData(char* subDataDesc, double** pBuffer, int* ht, int* wd, int band=1);

DLL_EXPORT int  ReadSubDataFloat(char* subDataDesc, float** pBuffer, int* ht, int* wd, int band=1);
DLL_EXPORT int  ReadSubDataShort(char* subDataDesc, short** pBuffer, int* ht, int* wd, int band=1);
DLL_EXPORT int  ReadSubDataUShort(char* subDataDesc, unsigned short** pBuffer, int* ht, int* wd, int band=1);

DLL_EXPORT int  ReadSubDataBandSum(char* subDataDesc);

DLL_EXPORT int  ReadSubDataCaliPara(char* subDataDesc, vector<double>& scale, vector<double>& offset);

/*************************   for HDF  ********************************/




#endif

