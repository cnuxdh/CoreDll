

#include "stdio.h"
#include "string.h"
#include "stdlib.h"

#include "wavelet.h"

//#include "Corelib/commonfile.h"
//#include "Corelib/image.h"


#define BYTE	unsigned char
#define UINT	unsigned int

/*
BYTE** LoadFromBmp(const char* filename, UINT &width, UINT &height)
{
	//handle file head
	long offset;
	int i;
	FILE *fp;
	fp=fopen(filename,"rb");
	if(!fp)
	{
		return NULL;
	}
	fseek(fp, 0x0a, SEEK_SET); // 0x0436 = 1078 = 54+256*4	
	fread(&offset,1,4,fp);
	fseek(fp,0x12,SEEK_SET);
	fread(&width,1,4,fp);
	fseek(fp,0x16,SEEK_SET);
	fread(&height,1,4,fp);
	fseek(fp,offset,SEEK_SET);
	//read the contents of the image
	BYTE** org_img = new BYTE*[height];
	for(i = 0; i < height; i++)
	{
		org_img[height - i - 1] = new BYTE[width];
		fread(org_img[height - i - 1], 1, width, fp);
	}
	//close the file
	fclose(fp);
	return org_img;
}

void SaveToBmp(const char* filename, const char* new_name, BYTE** rec_img, const UINT new_width, const UINT new_height)
{
	//handle file head
	long offset;
	int i;
	FILE *fp, *fpr;
	fp=fopen(filename,"rb");
	fpr=fopen(new_name,"wb");
	fseek(fp, 0x0a, SEEK_SET); // 0x0436 = 1078 = 54+256*4	
	fread(&offset,1,4,fp);
	rewind(fp);
	//copy file head
	for(i = 0; i < offset; i++)
	{
		fputc(fgetc(fp), fpr);
	}
	//new size
	fseek(fpr,0x12,SEEK_SET);
	fwrite(&new_width,1,4,fpr);
	fseek(fpr,0x16,SEEK_SET);
	fwrite(&new_height,1,4,fpr);
	fseek(fpr,offset,SEEK_SET);
	//fill the contents of the image
	for(i = 0; i < new_height; i++)
	{
		fwrite(rec_img[new_height - i - 1], 1, new_width, fpr);
	}
	//close the file
	fclose(fp);
	fclose(fpr);
}

*/
float** Select_Filter(const char* wavename, UINT &fil_len)
{
	float** filter_set;
	filter_set = new float*[4];
	int i;
	float* de_filter;
	//low pass decomposition filter
	if(!strcmp(wavename, "haar"))
	{
		fil_len = 2;
		float data[2] = {0.7071, 0.7071};
		de_filter = new float[fil_len]; 
		for(i = 0; i < fil_len; i++)
		{
			de_filter[i] = data[i];
		}
	}	
	else if(!strcmp(wavename, "db2"))
	{
		fil_len = 4;
		float data[4] = {-0.1294, 0.2241, 0.8365, 0.4830}; 
		de_filter = new float[fil_len]; 
		for(i = 0; i < fil_len; i++)
		{
			de_filter[i] = data[i];
		}
	}
	else if(!strcmp(wavename, "db3"))
	{
		fil_len = 6;
		float data[6] = {0.0352, -0.0854, -0.1350, 0.4599, 0.8069, 0.3327};
		de_filter = new float[fil_len]; 
		for(i = 0; i < fil_len; i++)
		{
			de_filter[i] = data[i];
		}
	}
	else if(!strcmp(wavename, "db4"))
	{
		fil_len = 8;
		float data[8] = {-0.0106, 0.0329, 0.0308, -0.1870, -0.0280, 0.6309, 0.7148, 0.2304};
		de_filter = new float[fil_len]; 
		for(i = 0; i < fil_len; i++)
		{
			de_filter[i] = data[i];
		}
	}
	else if(!strcmp(wavename, "db5"))
	{
		fil_len = 10;
		float data[10] = {0.0033, -0.0126, -0.0062, 0.0776, -0.0322, -0.2423, 0.1384, 0.7243, 0.6038, 0.1601};
		de_filter = new float[fil_len]; 
		for(i = 0; i < fil_len; i++)
		{
			de_filter[i] = data[i];
		}
	}
	else
	{
		return NULL;
	}
	//filter set 
	for(i = 0; i < 4; i++)
	{
		filter_set[i] = new float[fil_len];
	} 
	for(i = 0; i < fil_len; i++)
	{
		filter_set[0][i] = de_filter[i];
		filter_set[1][i] = de_filter[fil_len - i - 1] * (i % 2 ? 1 : -1);
		filter_set[2][i] = de_filter[fil_len - i - 1];
		filter_set[3][i] = de_filter[i] * (i % 2 ? -1 : 1);
	}
	return filter_set;
}

BYTE PixelClip(int pixel)
{
	if(pixel < 0)
	{
		pixel = 0;
	}
	else if(pixel > 255)
	{
		pixel = 255;
	}
	return(BYTE(pixel));
}

float* Filtering(const float* input, const float* filter, const UINT inlen, const UINT filen)
{
	int i, j;
	float* output = new float[inlen];
	for(i = 0; i < inlen; i++)
	{
		output[i] = 0;
		for(j = 0; j <= (i < (filen - 1) ? i : filen - 1); j++)
		{
			output[i] += input[i - j] * filter[j];
		}
	}
	return output;
}


float** DWT(BYTE** org_img, const float* ld_filter, const float* hd_filter,
			const UINT &width, const UINT &height, const UINT &fil_len, const UINT level)
{
	int i, j, k;
	UINT ex_len = fil_len;
	for(i = 0; i < level; i++)
	{
		ex_len *= 2; 
	}
	UINT ex_width  = width + ex_len;
	UINT ex_height = height + ex_len;

	//normalization
	float** std_img = new float*[ex_height];
	float** dec_img = new float*[ex_height];
	BYTE**  disp_img = new BYTE*[ex_height];
	for(i = 0; i < ex_height; i++)
	{
		std_img[i] = new float[ex_width];
		dec_img[i] = new float[ex_width];
		disp_img[i] = new BYTE[ex_width];
		for(j = 0; j < ex_width; j++)
		{
			if(i < height && j < width)
			{
				std_img[i][j] = float(org_img[i][j]) / 256;		
			}
			else
			{
				std_img[i][j] = 0;
			}
			dec_img[i][j] = std_img[i][j];
			disp_img[i][j] = 0;
		}
	}
	//decomposition
	float **dec_img_h, **dec_img_v;
	float *ltemp, *htemp;
	for(k = 0; k < level; k++)
	{
		dec_img_h = new float*[ex_height];
		dec_img_v = new float*[ex_height];
		for(i = 0; i < ex_height; i++)
		{
			dec_img_h[i] = new float[ex_width];
			dec_img_v[i] = new float[ex_width];
			for(j = 0; j < ex_width; j++)
			{
				dec_img_h[i][j] = 0;
				dec_img_v[i][j] = 0;
			}
		}

		//horizontal filtering
		for(i = 0; i < ex_height; i++)
		{
			ltemp = Filtering(dec_img[i], ld_filter, ex_width, fil_len); //low filter
			htemp = Filtering(dec_img[i], hd_filter, ex_width, fil_len); //high filter
			
			//downsampling
			for(j = 0; j < ex_width / 2; j++)
			{
				dec_img_h[i][j] = ltemp[2 * j]; 
				dec_img_h[i][j + ex_width / 2] = htemp[2 * j];
			}
		}

		/*
		//////////////////////////////////////////////////////////////////////////
		//added by xdh, to save the image for debug, 2013.7.4
		for(i = 0; i < ex_height; i++)
		{	
			for(j = 0; j < ex_width; j++)
			{
				disp_img[i][j] = PixelClip(int(256 * dec_img_h[i][j]));
			}
		}
		SaveToBmp(filename, "d:\\horizontalFilter.bmp", disp_img, ex_width, ex_height);
		//////////////////////////////////////////////////////////////////////////
		*/

		//vertical filtering
		float* tran_temp = new float[ex_height];
		for(j = 0; j < ex_width; j++)
		{
			for(i = 0; i < ex_height; i++)
			{
				tran_temp[i] = dec_img_h[i][j];
			}
			ltemp = Filtering(tran_temp, ld_filter, ex_height, fil_len);
			htemp = Filtering(tran_temp, hd_filter, ex_height, fil_len);
			
			//downsampling
			for(i = 0; i < ex_height / 2; i++)
			{
				dec_img_v[i][j] = ltemp[2 * i];
				dec_img_v[i + ex_height / 2][j] = htemp[2 * i];
			}
		}

		//copy to dec_img
		for(i = 0; i < ex_height; i++)
		{	
			for(j = 0; j < ex_width; j++)
			{
				dec_img[i][j] = dec_img_v[i][j];
			}
		}
		//zoom in
		ex_width /= 2;
		ex_height /= 2;
	}

	//save the decomposite image
	for(i = 0; i < ex_height; i++)
	{
		for(j = 0; j < ex_width; j++)
		{			
			disp_img[i][j] = PixelClip(int(256 * dec_img[i][j]));
		}
	}
	
	/*
	char* new_name = "dec.bmp";
	SaveToBmp(filename, new_name, disp_img, ex_width, ex_height);
	printf("\nThe decomposite image has been saved to dec.bmp\n");
	*/
	
	return dec_img;
}


/* interface written by xiedonghai, 2013.7.5
*/
void DWT(unsigned char* org_img, int height, int width, float** pWavelet, int* dstHt, int* dstWd,  
		 const float* ld_filter, const float* hd_filter, const UINT &fil_len, const UINT level )
{
	int i, j, k;
	UINT ex_len = fil_len;
	for(i = 0; i < level; i++)
	{
		ex_len *= 2; 
	}
	UINT ex_width  = width + ex_len;
	UINT ex_height = height + ex_len;

	*dstHt = ex_height;
	*dstWd = ex_width;
	(*pWavelet) = (float*)malloc(ex_height*ex_width*sizeof(float));

	//malloc memory
	float*  tran_temp = new float[ex_height];
	float** std_img = new float*[ex_height];
	float** dec_img = new float*[ex_height];
	BYTE**  disp_img = new BYTE*[ex_height];
	float** dec_img_h = new float*[ex_height];
	float** dec_img_v = new float*[ex_height];
	float *ltemp, *htemp;
	for(i = 0; i < ex_height; i++)
	{
		std_img[i]  = new float[ex_width];
		dec_img[i]  = new float[ex_width];
		disp_img[i] = new BYTE[ex_width];
		for(j = 0; j < ex_width; j++)
		{
			if(i < height && j < width)
			{
				std_img[i][j] = float(org_img[i*width+j]) / 256.0;		
			}
			else
			{
				std_img[i][j] = 0;
			}
			dec_img[i][j] = std_img[i][j];
			disp_img[i][j] = 0;
		}
		
		dec_img_h[i] = new float[ex_width];
		dec_img_v[i] = new float[ex_width];
		for(j = 0; j < ex_width; j++)
		{
			dec_img_h[i][j] = 0;
			dec_img_v[i][j] = 0;
		}
	}
	//////////////////////////////////////////////////////////////////////////

	//SaveBmp("d:\\src.bmp", org_img, height, width);

	//decomposition
	for(k = 0; k < level; k++)
	{
		//horizontal filtering
		for(i = 0; i < ex_height; i++)
		{
			ltemp = Filtering(dec_img[i], ld_filter, ex_width, fil_len); //low filter
			htemp = Filtering(dec_img[i], hd_filter, ex_width, fil_len); //high filter
			//downsampling
			for(j = 0; j < ex_width / 2; j++)
			{
				dec_img_h[i][j] = ltemp[2 * j]; 
				dec_img_h[i][j + ex_width / 2] = htemp[2 * j];
			}
			delete[] ltemp;
			delete[] htemp;
		}	

		/*
		//////////////////////////////////////////////////////////////////////////
		//added by xdh, to save the image for debug, 2013.7.5
		for(i = 0; i < ex_height; i++)
		{	
			for(j = 0; j < ex_width; j++)
			{
				disp_img[i][j] = PixelClip(int(256 * dec_img_h[i][j]));
			}
		}
		//SaveToBmp(filename, "d:\\horizontalFilter.bmp", disp_img, ex_width, ex_height);
		SaveBmp("d:\\hfilter.bmp", disp_img, ex_height, ex_width);
		//////////////////////////////////////////////////////////////////////////
		*/
		
		//vertical filtering
		for(j = 0; j < ex_width; j++)
		{
			for(i = 0; i < ex_height; i++)
			{
				tran_temp[i] = dec_img_h[i][j];
			}
			ltemp = Filtering(tran_temp, ld_filter, ex_height, fil_len);
			htemp = Filtering(tran_temp, hd_filter, ex_height, fil_len);
			//downsampling
			for(i = 0; i < ex_height / 2; i++)
			{
				dec_img_v[i][j] = ltemp[2 * i];
				dec_img_v[i + ex_height / 2][j] = htemp[2 * i];
			}
			delete[] ltemp;
			delete[] htemp;
		}
		//copy to dec_img
		for(i = 0; i < ex_height; i++)
		{	
			for(j = 0; j < ex_width; j++)
			{
				dec_img[i][j] = dec_img_v[i][j];
			}
		}

		//zoom in
		ex_width /= 2;
		ex_height /= 2;
	}


	ex_width  = width + ex_len;
	ex_height = height + ex_len;
	//save the decomposite image
	for(i = 0; i < ex_height; i++)
	{
		for(j = 0; j < ex_width; j++)
		{			
			disp_img[i][j] = PixelClip(int(256 * dec_img[i][j]));
			(*pWavelet)[i*ex_width+j] = dec_img[i][j];
		}
	}
	//SaveBmp("d:\\wavelet.bmp", disp_img, ex_height, ex_width);

	//release memory
	for(i = 0; i < ex_height; i++)
	{
		delete[] std_img[i];
		delete[] disp_img[i];
		delete[] dec_img[i];
		delete[] dec_img_h[i];
		delete[] dec_img_v[i];
	}
	delete[] tran_temp;
	//////////////////////////////////////////////////////////////////////////
}

BYTE** IDWT(float** dec_img, const float* lr_filter, const float* hr_filter, 
			const UINT &width, 	const UINT &height, const UINT &fil_len, const UINT level)
{
	int i, j, k;
	UINT ex_len = fil_len;
	for(i = 0; i < level; i++)
	{
		ex_len *= 2; 
	}
	UINT ex_width = width + ex_len;
	UINT ex_height = height + ex_len;
	float** std_img = new float*[ex_height];
	for(i = 0; i < ex_height; i++)
	{
		std_img[i] = new float[ex_width];
		for(j = 0; j < ex_width; j++)
		{
			std_img[i][j] = dec_img[i][j];
		}
	}
	BYTE** rec_img = new BYTE*[height];
	for(i = 0; i < height; i++)
	{
		rec_img[i] = new BYTE[width];
		for(j = 0; j < width; j++)
		{
			rec_img[i][j] = 0;
		}
	}
	//reconstruction
	float **rec_img_h, **rec_img_v;
	float *ltemp, *htemp, *tran_temp_l, *tran_temp_h, *temp_l, *temp_h;
	for(i = 1; i < level; i++)
	{
		ex_width /= 2;
		ex_height /= 2;
	}
	for(k = 0; k < level; k++)
	{	
		rec_img_h = new float*[ex_height];
		rec_img_v = new float*[ex_height];
		for(i = 0; i < ex_height; i++)
		{
			rec_img_h[i] = new float[ex_width];
			rec_img_v[i] = new float[ex_width];
			for(j = 0; j < ex_width; j++)
			{
				rec_img_h[i][j] = 0;
				rec_img_v[i][j] = 0;
			}
		}
		//vertical filtering
		tran_temp_l = new float[ex_height];
		tran_temp_h = new float[ex_height];
		for(j = 0; j < ex_width; j++)
		{
			for(i = 0; i < ex_height / 2; i++)
			{
				//upsampling
				tran_temp_l[2 * i] = std_img[i][j];
				tran_temp_h[2 * i] = std_img[i + ex_height / 2][j];
				tran_temp_l[2 * i + 1] = 0;
				tran_temp_h[2 * i + 1] = 0;
			}
			ltemp = Filtering(tran_temp_l, lr_filter, ex_height, fil_len);
			htemp = Filtering(tran_temp_h, hr_filter, ex_height, fil_len);
			for(i = 0; i < ex_height; i++)
			{
				rec_img_v[i][j] = ltemp[i] + htemp[i];
			}
		}
		//horizontal filtering
		temp_l = new float[ex_width];
		temp_h = new float[ex_width];
		for(i = 0; i < ex_height; i++)
		{
			for(j = 0; j < ex_width / 2; j++)
			{
				//upsampling
				temp_l[2 * j] = rec_img_v[i][j];
				temp_h[2 * j] = rec_img_v[i][j + ex_width / 2];
				temp_l[2 * j + 1] = 0;
				temp_h[2 * j + 1] = 0;
			}
			ltemp = Filtering(temp_l, lr_filter, ex_width, fil_len);
			htemp = Filtering(temp_h, hr_filter, ex_width, fil_len);
			for(j = 0; j < ex_width; j++)
			{
				rec_img_h[i][j] = ltemp[j] + htemp[j];
			}
		}
		//copy to std_img
		for(i = 0; i < ex_height - fil_len + 1; i++)
		{	
			for(j = 0; j < ex_width - fil_len + 1; j++)
			{
				std_img[i][j] = rec_img_h[i + fil_len - 1][j + fil_len - 1];
			}
		}
		//zoom out
		ex_width *= 2;
		ex_height *= 2;
	}
	//unnormalization
	for(i = 0; i < height; i++)
	{
		for(j = 0; j < width; j++)
		{			
			rec_img[i][j] = PixelClip(int(256 * std_img[i][j]));
		}
	}
	return rec_img;
}


void DWTProcess(unsigned char* pImage, int ht, int wd, float** pDWT, int* dstHt, int* dstWd, int level)
{
	//select wavelet filters 
	float** filter_set;
	UINT fil_len;
	filter_set = Select_Filter("db2", fil_len);
	if(filter_set)
	{
		//printf("\nThe wavelet filter has been found. Length = %d\n", fil_len);
	}
	else
	{
		printf("\nSorry. The wavelet filter has not been found.\n");
		return;
	}
	float *ld_filter, *hd_filter, *lr_filter, *hr_filter;
	ld_filter = filter_set[0];  //low filter for decomposition
	hd_filter = filter_set[1];  //high filter for decompostion
	lr_filter = filter_set[2];
	hr_filter = filter_set[3];

	//out filter
	printf("Low filter: \n");
	for(int i=0; i<fil_len; i++)
	{
		printf("%lf ", ld_filter[i]);
	}
	printf("\n");
	printf("High filter: \n");
	for(int i=0; i<fil_len; i++)
	{
		printf("%lf ", hd_filter[i]);
	}
	printf("\n");

	//DFT
	DWT(pImage, ht, wd, pDWT, dstHt, dstWd, ld_filter, hd_filter, fil_len, level);
}

void main_test()
{
	/*
	//load the image file
	char* filename = new char[20];
	printf("\nPlease enter the image filename(bmp):");
	gets(filename);
	if(!org_img)
	{
		printf("\nSorry. The file does not exist.\n");
		return;
	}
	else
	{
		printf("\nA %d*%d file has been loaded.\n", height, width);
	}	
	//select the wavelet filter
	char* wavename = new char[20];
	printf("\nPlease enter the name of wavelet filter:");
	gets(wavename);
	*/

	char* filename = "d:\\data\\lena-gray.bmp";
	char* wavename = "db2"; 
    printf("%s \n", filename);

	BYTE** org_img;
	UINT width, height;
	//org_img = LoadFromBmp(filename, width, height);
	printf("%d %d \n", height, width);

		
	float** filter_set;
	UINT fil_len;
	filter_set = Select_Filter(wavename, fil_len);
	if(filter_set)
	{
		printf("\nThe wavelet filter has been found. Length = %d\n", fil_len);
	}
	else
	{
		printf("\nSorry. The wavelet filter has not been found.\n");
		return;
	}
	float *ld_filter, *hd_filter, *lr_filter, *hr_filter;
	ld_filter = filter_set[0];  //low filter for decomposition
	hd_filter = filter_set[1];  //high filter for decompostion
	lr_filter = filter_set[2];
	hr_filter = filter_set[3];
	
	//wavelet decomposition
	UINT level;
	printf("\nPlease enter the transform level:");
	scanf("%d", &level);

	float** dec_img;
	dec_img = DWT(org_img, ld_filter, hd_filter, width, height, fil_len, level);	
	//wavelet reconstruction
	BYTE** rec_img;
	rec_img = IDWT(dec_img, lr_filter, hr_filter, width, height, fil_len, level);
	
	//save the image file
	//char* new_name = "rec.bmp";
	//SaveToBmp(filename, new_name, rec_img, width, height);
	//printf("\nThe reconstruct image has been saved to rec.bmp\n");
}
