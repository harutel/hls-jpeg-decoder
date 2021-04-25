/***************************************************************************/
/*                                                                         */
/*  File: main.cpp                                                         */
/*  Autor: bkenwright@xbdev.net                                            */
/*  URL: www.xbdev.net                                                     */
/*                                                                         */
/***************************************************************************/
/*
	Jpeg File Format Explained
*/
/***************************************************************************/

//#include <windows.h>

#include <stdio.h>		// sprintf(..), fopen(..)
#include <stdarg.h>     // So we can use ... (in dprintf)
#include <string.h>		// memset(..)

#include "openjpg.h"	// OpenJpgFile(..)
#include "loadjpg.h"	// ConvertJpgFile(..)

//#include "savejpg.h"    // SaveJpgFile(..)
//#include "closejpg.h"    // SaveJpgFile(..)

/***************************************************************************/
/*                                                                         */
/* FeedBack Data                                                           */
/*                                                                         */
/***************************************************************************/

//Saving debug information to a log file
void dprintf(const char *fmt, ...) 
{
/*
	va_list parms;
	char buf[256];

	// Try to print in the allocated space.
	va_start(parms, fmt);
	vsprintf (buf, fmt, parms);
	va_end(parms);

	// Write the information out to a txt file
	FILE *fp = fopen("output.txt", "a+");
	fprintf(fp, "%s", buf);
	fclose(fp);
*/
}// End dprintf(..)


 
/***************************************************************************/
/*                                                                         */
/* Entry Point                                                             */
/*                                                                         */
/***************************************************************************/


//int __stdcall WinMain (HINSTANCE hInst, HINSTANCE hPrev, LPSTR lpCmd, int nShow)
//{
//
//	// Create a jpg from a bmp
//	//SaveJpgFile("smiley.bmp", "ex.jpg");
//
//	// Create a bmp from a jpg
//	ConvertJpgFile("smiley.jpg", "smiley.bmp");
//
//	//ConvertJpgFile("imageJPEG3_1.jpg", "imageJPEG3_1.bmp");
//	//ConvertJpgFile("imageJPEG3_2.jpg", "imageJPEG3_2.bmp");
//
//	//ConvertJpgFile("test.jpg", "test.bmp");
//
//	//ConvertJpgFile("cross.jpg", "cross.bmp");
//
//	return 0;
//}// End WinMain(..)

/***************************************************************************/
//
// Take Jpg data, i.e. jpg file read into memory, and decompress it to an
// array of rgb pixel values.
//
// Note - Memory is allocated for this function, so delete it when finished
//
///***************************************************************************/
//int DecodeJpgFileData(const unsigned char* buf, // Jpg file in memory
//					  const int sizeBuf,		// Size jpg in bytes in memory
//					  unsigned char** rgbpix,	// Output rgb pixels
//					  unsigned int* width,		// Output image width
//					  unsigned int* height)		// Output image height
//{
//	// Allocate memory for our decoded jpg structure, all our data will be
//	// decompressed and stored in here for the various stages of our jpeg decoding
//	stJpegData* jdec = new stJpegData();
//	memset(jdec, 0, sizeof(stJpegData));
//	 printf("size:%d",sizeof(stJpegData));
//	if (jdec == NULL)
//	{
//		dprintf("Not enough memory to alloc the structure need for decompressing\n");
//		///DBG_HALT;
//		return 0;
//	}
//
//	// Start Parsing.....reading & storing data
//	if (JpegParseHeader(jdec, buf, sizeBuf)<0)
//	{
//		dprintf("ERROR > parsing jpg header\n");
//		///DBG_HALT;
//	}
//
//	// We've read it all in, now start using it, to decompress and create rgb values
//	dprintf("Decoding JPEG image...\n");
//	JpegDecode(jdec);
//
//	// Get the size of the image
//	JpegGetImageSize(jdec, width, height);
//
//	*rgbpix = jdec->m_rgb;
//
//	// Release the memory for our jpeg decoder structure jdec
//	delete jdec;
//
//	return 1;
//}

///***************************************************************************/
////
////  Returns the size of the file in bytes
////
///***************************************************************************/
//inline int FileSize(FILE *fp)
//{
//	long pos;
//	fseek(fp, 0, SEEK_END);
//	pos = ftell(fp);
//	fseek(fp, 0, SEEK_SET);
//	return pos;
//}

///***************************************************************************/
//int DecodeJpgFileData(const unsigned char* buf, // Jpg file in memory
//					  const int sizeBuf,		// Size jpg in bytes in memory
//					  unsigned char** rgbpix,	// Output rgb pixels
//					  unsigned int* width,		// Output image width
//					  unsigned int* height)		// Output image height
//{
//	// Allocate memory for our decoded jpg structure, all our data will be
//	// decompressed and stored in here for the various stages of our jpeg decoding
//	stJpegData* jdec = new stJpegData();
//	memset(jdec, 0, sizeof(stJpegData));
//	 printf("size:%d",sizeof(stJpegData));
//	if (jdec == NULL)
//	{
//		dprintf("Not enough memory to alloc the structure need for decompressing\n");
//		///DBG_HALT;
//		return 0;
//	}
//
//	// Start Parsing.....reading & storing data
//	if (JpegParseHeader(jdec, buf, sizeBuf)<0)
//	{
//		dprintf("ERROR > parsing jpg header\n");
//		///DBG_HALT;
//	}
//
//	// We've read it all in, now start using it, to decompress and create rgb values
//	dprintf("Decoding JPEG image...\n");
//	JpegDecode(jdec);
//
//	// Get the size of the image
//	JpegGetImageSize(jdec, width, height);
//
//	*rgbpix = jdec->m_rgb;
//
//	// Release the memory for our jpeg decoder structure jdec
//	delete jdec;
//
//	return 1;
//}

///***************************************************************************/
////
//// Save a buffer in 24bits Bitmap (.bmp) format
////
///***************************************************************************/
//inline void WriteBMP24(const char* szBmpFileName, int Width, int Height, unsigned char* RGB)
//{
//	#pragma pack(1)
//	struct stBMFH // BitmapFileHeader & BitmapInfoHeader
//	{
//		// BitmapFileHeader
//		char         bmtype[2];     // 2 bytes - 'B' 'M'
//		unsigned int iFileSize;     // 4 bytes
//		short int    reserved1;     // 2 bytes
//		short int    reserved2;     // 2 bytes
//		unsigned int iOffsetBits;   // 4 bytes
//		// End of stBMFH structure - size of 14 bytes
//		// BitmapInfoHeader
//		unsigned int iSizeHeader;    // 4 bytes - 40
//		unsigned int iWidth;         // 4 bytes
//		unsigned int iHeight;        // 4 bytes
//		short int    iPlanes;        // 2 bytes
//		short int    iBitCount;      // 2 bytes
//		unsigned int Compression;    // 4 bytes
//		unsigned int iSizeImage;     // 4 bytes
//		unsigned int iXPelsPerMeter; // 4 bytes
//		unsigned int iYPelsPerMeter; // 4 bytes
//		unsigned int iClrUsed;       // 4 bytes
//		unsigned int iClrImportant;  // 4 bytes
//		// End of stBMIF structure - size 40 bytes
//		// Total size - 54 bytes
//	};
//	#pragma pack()
//
//	// Round up the width to the nearest DWORD boundary
//	int iNumPaddedBytes = 4 - (Width * 3) % 4;
//	iNumPaddedBytes = iNumPaddedBytes % 4;
//
//	stBMFH bh;
//	memset(&bh, 0, sizeof(bh));
//	bh.bmtype[0]='B';
//	bh.bmtype[1]='M';
//	bh.iFileSize = (Width*Height*3) + (Height*iNumPaddedBytes) + sizeof(bh);
//	bh.iOffsetBits = sizeof(stBMFH);
//	bh.iSizeHeader = 40;
//	bh.iPlanes = 1;
//	bh.iWidth = Width;
//	bh.iHeight = Height;
//	bh.iBitCount = 24;
//
//
//	char temp[2048]={0};
//	sprintf(temp, "%s", szBmpFileName);
//	FILE* fp = fopen(temp, "wb");
//	///DBG_ASSERT( fp ); // Error creating file - valid file name?  valid location?
//	fwrite(&bh, sizeof(bh), 1, fp);
//	for (int y=Height-1; y>=0; y--)
//	{
//		for (int x=0; x<Width; x++)
//		{
//			int i = (x + (Width)*y) * 3;
//			unsigned int rgbpix = (RGB[i]<<16)|(RGB[i+1]<<8)|(RGB[i+2]<<0);
//			fwrite(&rgbpix, 3, 1, fp);
//		}
//		if (iNumPaddedBytes>0)
//		{
//			unsigned char pad = 0;
//			fwrite(&pad, iNumPaddedBytes, 1, fp);
//		}
//	}
//	fclose(fp);
//}
//
///***************************************************************************/
////
//// Load one jpeg image, and decompress it, and save the result.
////
///***************************************************************************/
//int ConvertJpgFile(char* szJpgFileInName, char * szBmpFileOutName)
//{
//	FILE *fp;
//	unsigned int lengthOfFile;
//	unsigned char *buf;
//
//
//	// Load the Jpeg into memory
//	fp = fopen(szJpgFileInName, "rb");
//	if (fp == NULL)
//	{
//		dprintf("Cannot open jpg file: %s\n", szJpgFileInName);
//		///DBG_HALT;
//		return 0;
//	}
//
//	lengthOfFile = FileSize(fp);
//	buf = new unsigned char[lengthOfFile + 4];// +4 is safety padding
//	if (buf == NULL)
//	{
//		dprintf("Not enough memory for loading file\n");
//		///DBG_HALT;
//		return 0;
//	}
//	fread(buf, lengthOfFile, 1, fp);
//	fclose(fp);
//
//	unsigned char* rgbpix = NULL;
//	unsigned int width  = 0;
//	unsigned int height = 0;
//	DecodeJpgFileData(buf, lengthOfFile, &rgbpix, &width, &height);
//
//	if (rgbpix==NULL)
//	{
//		dprintf("Failed to decode jpg\n");
//		///DBG_HALT;
//		return 0;
//	}
//
//	// Delete our data we read in from the file
//	delete[] buf;
//
//	// Save it
//	WriteBMP24(szBmpFileOutName, width, height, rgbpix);
//
//	// Since we don't need the pixel information anymore, we must
//	// release this as well
//	delete[] rgbpix;
//
//	return 1;
//}

int main()
{
	char jpg_in[] = "/mnt/FILES/p/hls-jpeg-decoder/data/Lenna.jpg";
	char bmp_out[] = "/mnt/FILES/p/hls-jpeg-decoder/out.bmp";

	ConvertJpgFile(jpg_in, bmp_out);

	return (0);
}


