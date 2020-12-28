// Sinc2D.cpp : Defines the entry point for the console application.
//

// Sinc2D (SincLin2Resize trapezoid kernel weighted) upsizer to iMul ratio.
// interpolate in-between input samples and re-assign to new integer size

#include "stdafx.h"
#include "windows.h"
#include "commctrl.h"

#include "tiffio.h"


#include <cmath>
#include "stdlib.h"


BYTE *g_pImageBuffer = 0, *g_pElImageBuffer, *g_pFilteredImageBuffer = 0;

float *g_pfImageBuffer = 0, *g_pfFilteredImageBuffer = 0;

float *g_pfKernel = 0;

int iKernelSize;
int iMul = 4;
int iTaps = 15;

int iWidth;
int iHeight;

//-----------------------------------------------------------------------------
// Msg: Display an error if needed
//-----------------------------------------------------------------------------
void Msg(TCHAR *szFormat, ...)
{
    printf(szFormat);
}

BOOL LoadTIFF8(void)
{
	TIFF *in_tiff;
	unsigned short bps, spp, cmpr;
	tstrip_t strip;
	uint32* bc;
	uint32 stripsize, rowspstrip;
	uint32 uiWidth, uiHeight;
	uint32 row, col;

	unsigned char *ucSrc_data;
	unsigned char *ucImageBuffer = (unsigned char*)g_pImageBuffer;

	// Open the src TIFF image
	if((in_tiff = TIFFOpen("in.tif", "r")) == NULL)
	{
		Msg(TEXT("Could not open src TIFF image"));
		return FALSE;
	}
	
	// Check that it is of a type that we support
	if((TIFFGetField(in_tiff, TIFFTAG_COMPRESSION, &cmpr) ==0 ) || (cmpr != COMPRESSION_NONE))
	{
		Msg(TEXT("Compressed TIFFs not supported"));
		TIFFClose(in_tiff);
		return FALSE;
	}

	if((TIFFGetField(in_tiff, TIFFTAG_BITSPERSAMPLE, &bps) == 0) || (bps != 8))
	{
		Msg(TEXT("Either undefined or unsupported number of bits per sample - must be 8"));
		TIFFClose(in_tiff);
		return FALSE;
	}

	if((TIFFGetField(in_tiff, TIFFTAG_SAMPLESPERPIXEL, &spp) == 0) || (spp != 3))
	{
		Msg(TEXT("Either undefined or unsupported number of samples per pixel - must be 3"));
		TIFFClose(in_tiff);
		return FALSE;
	}

	TIFFGetField(in_tiff, TIFFTAG_IMAGEWIDTH, &uiWidth);
	iWidth = uiWidth;

	TIFFGetField(in_tiff, TIFFTAG_IMAGELENGTH, &uiHeight);
	iHeight = uiHeight;

	// finally allocate mem bufs 
	g_pImageBuffer = (BYTE*)malloc(iWidth * iHeight * 3);
	g_pFilteredImageBuffer = (BYTE*)malloc(iWidth * iHeight * iMul * iMul * 3);

	ZeroMemory(g_pFilteredImageBuffer, (iWidth * iHeight * iMul * iMul * 3));

	g_pfImageBuffer = (float*)malloc(iWidth * iHeight * 3 * sizeof(float));
	ZeroMemory(g_pImageBuffer, (iWidth * iHeight * 3));


	TIFFGetField(in_tiff, TIFFTAG_STRIPBYTECOUNTS, &bc);
	stripsize = bc[0];
	// test if all strips are equal
	uint32 allstripsize = 0;

	for (strip = 0; strip < TIFFNumberOfStrips(in_tiff); strip++) 
	{
		if (bc[strip] != stripsize) 
		{
			Msg(TEXT("Strip sizes unequal"));
//			return FALSE;
		}

		allstripsize += bc[strip];
	}

	ucSrc_data = (unsigned char*)malloc(allstripsize);

	TIFFGetField(in_tiff, TIFFTAG_ROWSPERSTRIP, &rowspstrip);

	// load from tiff to temp buffer
	int iBytesRead = 0;
	for(row = 0; row < uiHeight; row+=rowspstrip)
	{

		uint32 readsize = stripsize;

		if (row*uiWidth * 3 + stripsize > allstripsize) continue; // temp workaround, need allocate sumof(stripsize) mem

		iBytesRead += TIFFReadRawStrip(in_tiff, row/rowspstrip, &ucSrc_data[row*uiWidth*3], stripsize);
	}

	printf("Read bytes from tiff %d , all image bytes must be %d\r\n", iBytesRead, iHeight*iWidth*3);

	// separate planes
	for(row=0; row < uiHeight; row++)
	{
		for(col=0; col < uiWidth; col++)
		{
			g_pImageBuffer[(row*uiWidth + col)] = ucSrc_data[(row*uiWidth+col)*3+0]; // R
			g_pImageBuffer[(uiWidth*uiHeight + row*uiWidth + col)] = ucSrc_data[(row*uiWidth+col)*3+1]; // G
			g_pImageBuffer[(uiWidth*uiHeight*2+ row*uiWidth + col)] = ucSrc_data[(row*uiWidth+col)*3+2]; // B
		}
	}

	free(ucSrc_data);
	TIFFClose(in_tiff);

	return TRUE;
}

int SaveTIFF8()
{
	TIFF *out_tiff;
	// Open the dst TIFF image
	if((out_tiff = TIFFOpen("out.tif", "w")) == NULL)
	{
		fprintf(stderr, "Could not open dst TIFF image\n");
		return(1);
	}


	TIFFSetField(out_tiff, TIFFTAG_IMAGEWIDTH, iWidth*iMul);
	TIFFSetField(out_tiff, TIFFTAG_BITSPERSAMPLE, 8);
	TIFFSetField(out_tiff, TIFFTAG_SAMPLESPERPIXEL, 3);
	TIFFSetField(out_tiff, TIFFTAG_ROWSPERSTRIP, 1);
	TIFFSetField(out_tiff, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
	TIFFSetField(out_tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
	TIFFSetField(out_tiff, TIFFTAG_FILLORDER, FILLORDER_MSB2LSB);
	TIFFSetField(out_tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);

	TIFFSetField(out_tiff, TIFFTAG_XRESOLUTION, 300);
	TIFFSetField(out_tiff, TIFFTAG_YRESOLUTION, 300);
	TIFFSetField(out_tiff, TIFFTAG_RESOLUTIONUNIT, RESUNIT_INCH);

	TIFFSetField(out_tiff, TIFFTAG_IMAGELENGTH, iHeight*iMul);

	unsigned char *ucDst_data  = (unsigned char*)malloc(iWidth*3*iMul);

	for (int row=0; row < iHeight*iMul; row++) // lines counter
		{

			for (int i = 0; i < iWidth*iMul; i++)
			{
				ucDst_data[i*3] = g_pFilteredImageBuffer[(row*iWidth*iMul + i)]; // R channel
				ucDst_data[i*3+1] = g_pFilteredImageBuffer[(iWidth*iHeight * iMul*iMul + row * iWidth * iMul + i)]; // G channel
				ucDst_data[i * 3 + 2] = g_pFilteredImageBuffer[(iWidth*iHeight * 2 * iMul*iMul + row*iWidth*iMul + i)]; // B channel
			}
		
			TIFFWriteRawStrip(out_tiff, row, ucDst_data, /*g_stripsize*/iWidth*3*iMul);
	}

	TIFFClose(out_tiff);

	free (ucDst_data);

	return 0;

}


void fill2DKernel(void)
{
	float fPi = 3.14159265358979f;
	int i,j;

	// our 2D just iMul-finely sampled kernel defined in output size

	/* Sinc (SincLin2) kernel 
	for (i = 0; i < iKernelSize; i++)
	{
		for (j = 0; j < iKernelSize; j++)
		{
			float fDist = sqrtf((float(iKernelSize / 2) - j)*(float(iKernelSize / 2) - j) + (float(iKernelSize / 2) - i)*(float(iKernelSize / 2) - i));

			// make kernel round in 2d
			if (fDist > iKernelSize / 2) continue;

			float fArg = (fPi*fDist/iMul); 

			if (fDist <= (iKernelSize / 4))
			{
				if (fArg != 0)
				{
					g_pfKernel[i*iKernelSize + j] = sinf(fArg) / fArg;
				}
				else
					g_pfKernel[i*iKernelSize + j] = 1.0f;

			}
			else if (fDist > (iKernelSize / 4) && fDist <= (iKernelSize/2))
			{
				g_pfKernel[i*iKernelSize + j] = (sinf(fArg) / fArg) * ((2 - (4 * fDist / iKernelSize))); // SincLin2Resize trapezoidal weighting
			}
			else
				g_pfKernel[i*iKernelSize + j] = 0.0f;
		} //j
	} //i
	*/

	/* Jinc weighted by Jinc - EWA Lanczos kernel */
	for (i = 0; i < iKernelSize; i++)
	{

		if (i == (iKernelSize / 2))
		{
			int dbr = 0;
		}
		for (j = 0; j < iKernelSize; j++)
		{
			float fDist = sqrtf((float(iKernelSize / 2) - j) * (float(iKernelSize / 2) - j) + (float(iKernelSize / 2) - i) * (float(iKernelSize / 2) - i));

			// make kernel round in 2d
			if (fDist > iKernelSize / 2) continue;

			float fArg = ( fPi * fDist / iMul);

			float fArg_w = fDist * 3.9f / (iKernelSize / 2);

			if (fArg != 0)
			{
				float fBess = 2.0f * std::cyl_bessel_jf(1, fArg) / fArg;
				float fW; // Jinc window
				if (fArg_w != 0)
				{
					fW = 2.0f * std::cyl_bessel_jf(1, fArg_w) / fArg_w;
				}
				else
				{
					fW = 1.0f;
				}
				
				g_pfKernel[i * iKernelSize + j] = fBess * fW;
			}
			else
				g_pfKernel[i * iKernelSize + j] = 1.0f;


		} //j
	} //i
	
	
	// normalize to 1
	float fSum = 0.0f;
	for (i = 0; i < iKernelSize; i++)
	{
		for (j = 0; j < iKernelSize; j++)
		{
			fSum += g_pfKernel[i*iKernelSize + j];
		}
	}

	for (i = 0; i < iKernelSize; i++)
	{
		for (j = 0; j < iKernelSize; j++)
		{
			g_pfKernel[i*iKernelSize + j] /= fSum;
			g_pfKernel[i*iKernelSize + j] *= (iMul * iMul); // energy dissipated at iMul^2 output samples, so 1 norm * iMul^2
		}
	}

	// fill 256 kernel images weighted by 8bit unsigned int
	for (int iSample = 2; iSample < 256; iSample++)
	{
		float* pfCurrKernel = g_pfKernel + iSample * iKernelSize * iKernelSize;
		// start from 2, 1st is kernel * 1
		for (i = 0; i < iKernelSize; i++)
		{
			for (j = 0; j < iKernelSize; j++)
			{
				pfCurrKernel[i * iKernelSize + j] = g_pfKernel[i * iKernelSize + j] * iSample;				
			}
		}
	}

	// TO DO : 1. Make kernel half height. 2. Add run-length encoding for non-zero line length.

}

void KernelProc(void)
{
	int row,col, k_row,k_col;

	int iWidthEl = iWidth + 2 * iKernelSize, iHeightEl = iHeight + 2 * iKernelSize;

	g_pElImageBuffer = (BYTE*)malloc(iWidthEl * iHeightEl);
	ZeroMemory(g_pElImageBuffer, iWidthEl * iHeightEl);

	g_pfFilteredImageBuffer = (float*)malloc(iWidthEl * iHeightEl * iMul * iMul * sizeof(float));
	ZeroMemory(g_pfFilteredImageBuffer, (iWidthEl * iHeightEl * iMul * iMul * sizeof(float)));

	// R-channel

	// fill center
	for (row=0; row < iHeight; row++) 
	{
		for (col = 0; col < iWidth; col++)
		{
			g_pElImageBuffer[(row + iKernelSize) * iWidthEl + (col + iKernelSize)] = g_pImageBuffer[(row*iWidth + col)]; // R channel
		}		
	}
	
	// fill upper strip - dup 1st line
	for (row = 0; row < iKernelSize; row++) 
	{
		for (col = 0; col < iWidth; col++)
		{
			g_pElImageBuffer[row * iWidthEl + (col + iKernelSize)] = g_pImageBuffer[col]; // R channel
		}
	}
	
	/* somwhere bug here in to wr buffer
	// fill lower strip - dup last line
	for (row = iHeight; row < iHeightEl; row++) 
	{
		for (col = 0; col < iWidth; col++)
		{
			g_pElImageBuffer[(row + iKernelSize) * iWidthEl + (col + iKernelSize)] = g_pImageBuffer[((iHeight - 1) * iWidth + col)]; // R channel
		}
	}
	*/
//	printf("Just %d row of %d starting R-channel. Be patient...\r\n", row, iHeightEl);
	
	// fill left strip - dup first column
	for (row = 0; row < iHeight; row++) 
	{
		for (col = 0; col < iKernelSize; col++)
		{
			g_pElImageBuffer[(row + iKernelSize) * iWidthEl + col] = g_pImageBuffer[(row * iWidth) + 0]; // R channel
		}
	}

	// fill right strip - dup last column
	for (row = 0; row < iHeight; row++) 
	{
		for (col = iWidth; col < (iWidth + iKernelSize); col++)
		{
			g_pElImageBuffer[(row + iKernelSize) * iWidthEl + col + iKernelSize] = g_pImageBuffer[(row * iWidth) + (iWidth - 1)]; // R channel
		}
	}
	

	/*
	// 2d convolution pass - mul to kernel
	for (row=0+iKernelSize/2; row < iHeight-iKernelSize/2; row++) // input lines counter
	{
		for (col = 0+iKernelSize/2; col < iWidth-iKernelSize/2; col++) // input cols counter
		{
			float fInpSample=g_pfImageBuffer[(row*iWidth + col)];

			for (k_row = 0; k_row < iKernelSize; k_row++)
			{
				for (k_col = 0; k_col < iKernelSize; k_col++)
				{
					float fKrn = g_pfKernel[(k_row)*iKernelSize + k_col];
					
					g_pfFilteredImageBuffer[(row*iMul - (iKernelSize / 2) + k_row)*iWidth*iMul + (col - (iKernelSize / 2))*iMul + k_col] += fInpSample*fKrn;
				} // k_col
			} // k_row

		} // col
	//	printf("Just %d row of %d done R-channel. Be patient...\r",row,iHeight);
	} // row
	*/
	
	// 2d convolution pass - read kernel lut and add row to row
	// process center 
#pragma omp parallel for
	for (row=0+(iTaps); row < iHeightEl-(iTaps); row++) // input lines counter
	{
		for (col = 0+(iTaps); col < iWidthEl-(iTaps); col++) // input cols counter
		{
			UINT uiInpSample=g_pElImageBuffer[(row*iWidthEl + col)]; // R-channel

			int iOutWidth = iWidthEl * iMul;

			// fast skip zero
			if (uiInpSample == 0) continue;

			float* pfCurrKernel = g_pfKernel + ((uiInpSample - 1) * iKernelSize * iKernelSize);

			float* pfFilteredImageBuffer = g_pfFilteredImageBuffer + ((row * iMul - (iKernelSize / 2)) * iOutWidth) + (col * iMul - (iKernelSize / 2));

			for (k_row = 0; k_row < iKernelSize; k_row++)
			{
				// add full kernel row to output
				for (k_col = 0; k_col < iKernelSize; k_col++)
				{
					pfFilteredImageBuffer[k_col] += pfCurrKernel[k_col];
				} // k_col
				pfFilteredImageBuffer += iOutWidth; // point to next start point in output buffer now
				pfCurrKernel += iKernelSize; // point to next kernel row now
			} // k_row

		} // col
	} // row
		
    // R-channel
	for (row = 0 ; row < iHeight * iMul; row++) 
	{
		for (col = 0; col < iWidth * iMul; col++)
		{
			unsigned char ucVal;
			float fVal = g_pfFilteredImageBuffer[(row + iKernelSize*iMul)*iWidthEl*iMul + col + iKernelSize * iMul]; // R channel

			fVal += 0.5f;

			if (fVal > 255.0f) 
			{
				fVal= 255.0f;
			}
			if (fVal < 0.0f) 
			{
				fVal = 0.0f;
			}
			ucVal = (unsigned char) fVal;

			g_pFilteredImageBuffer[(row*iWidth*iMul + col)] = ucVal; // R channel
			g_pFilteredImageBuffer[(iWidth*iHeight*iMul*iMul + row*iWidth*iMul + col)] = ucVal; //temp to G 
			g_pFilteredImageBuffer[(iWidth*iHeight*2*iMul*iMul + row*iWidth*iMul + col)] = ucVal; //temp to B 
		}		
	}
	

}

int main(int argc, char* argv[])
{
	printf("Hello World! First param is mul, second taps.\n");


	if (argc > 1)
	{
		argc--;
		argv++;
		iMul = atoi(argv[0]);

		if (argc > 1)
		{
			argc--;
			argv++;
			iTaps = atoi(argv[0]);
		}
	}

	printf("in.tif mul=%d taps=%d\n",iMul,iTaps);

	iKernelSize = (iMul*iTaps*2);

	// size of 2d kernel array iMul*iMul of KernelSize*KernelSize
	// add kernel LUT for 8bit unsigned input
	g_pfKernel = (float*)malloc(iKernelSize * iKernelSize * sizeof(float) * 256);
	

	ZeroMemory(g_pfKernel, iKernelSize * iKernelSize * sizeof(float) * 256);

	fill2DKernel();

	if (LoadTIFF8() == FALSE)
		goto end;

	printf("Processing input w=%d h=%d\r\n", iWidth, iHeight);

//	KernelProc();

	// test perf
// Measure Performance
	printf("\n\nPerformance Measurement: \n");
	int numIter = 5, iter;
	LARGE_INTEGER proc_freq;
	QueryPerformanceFrequency(&proc_freq);
	LARGE_INTEGER startTime;
	QueryPerformanceCounter(&startTime);

	float times[20000];
	__int64 starttime;
	__int64 endtime;

	for (iter = 0; iter < numIter; iter++)
	{
		QueryPerformanceCounter((LARGE_INTEGER*)&starttime);

		KernelProc();

		QueryPerformanceCounter((LARGE_INTEGER*)&endtime);
		endtime -= starttime;
		times[iter] = (float)endtime / (float)proc_freq.QuadPart;
	}

	float minval = 100000000;
	float maxval = 0;
	float average = 0;
	for (iter = 0; iter < numIter; iter++)
	{
		if (times[iter] > maxval) maxval = times[iter];
		if (times[iter] < minval) minval = times[iter];

		average += times[iter];
	}

	average /= numIter;

	printf("\n min: %g max: %g avg: %g (msec?)\n", minval * 1000, maxval * 1000, average * 1000);
	
	

	SaveTIFF8();
end:
	free (g_pImageBuffer);
	free(g_pElImageBuffer);
    free (g_pFilteredImageBuffer);

	free (g_pfImageBuffer);
	free (g_pfFilteredImageBuffer);

	free (g_pfKernel);

	return 0;
}

