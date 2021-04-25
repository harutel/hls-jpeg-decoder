	/***************************************************************************/
/*                                                                         */
/*  File: loadjpg.cpp                                                      */
/*  Author: bkenwright@xbdev.net                                           */
/*  Date: 19-01-06                                                         */
/*                                                                         */
/*  Revised: 26-07-07                                                      */
/*                                                                         */
/***************************************************************************/
/*
    About:
	Simplified jpg/jpeg decoder image loader - so we can take a .jpg file
	either from memory or file, and convert it either to a .bmp or directly
	to its rgb pixel data information.

	Simplified, and only deals with basic jpgs, but it covers all the
	information of how the jpg format works :)

	Can be used to convert a jpg in memory to rgb pixels in memory.

	Or you can pass it a jpg file name and an output bmp filename, and it
	loads and writes out a bmp file.

	i.e.
	ConvertJpgFile("cross.jpg", "cross.bmp")
*/
/***************************************************************************/

//#pragma once

#include <stdio.h>		// sprintf(..), fopen(..)
#include <stdarg.h>     // So we can use ... (in dprintf)
#include <string.h>		// memset(..)
#include <math.h>		// sqrt(..), cos(..)
#include <stdint.h>

#include "loadjpg.h"
#include "openjpg.h"

//#include "hls_stream.h"	///


//int ll=0;

#define JPG_SPEEDUP

extern void dprintf(const char *fmt, ...);

/*
__forceinline void dprintf(const char *fmt, ...)
{
	va_list parms;
	char buf[256];

	// Try to print in the allocated space.
	va_start(parms, fmt);
	vsprintf (buf, fmt, parms);
	va_end(parms);

	// Write the information out to a txt file
	FILE *fp = fopen("output.txt", "a+");
	fprintf(fp, "%s",  buf);
	fclose(fp);

}// End dprintf(..)
*/

/***************************************************************************/

//#define DQT 	 0xDB	// Define Quantization Table
//#define SOF 	 0xC0	// Start of Frame (size information)
//#define DHT 	 0xC4	// Huffman Table
//#define SOI 	 0xD8	// Start of Image
//#define SOS 	 0xDA	// Start of Scan
//#define EOI 	 0xD9	// End of Image, or End of File
//#define APP0	 0xE0
//
//#define BYTE_TO_WORD(x) (((x)[0]<<8)|(x)[1])


//const int HUFFMAN_TABLES	=	4;
//const int COMPONENTS		=	4;

//const int  cY	 = 1;
//const int  cCb	 = 2;
//const int  cCr	 = 3;


static int ZigZagArray[64] =
{
	0,   1,   5,  6,   14,  15,  27,  28,
	2,   4,   7,  13,  16,  26,  29,  42,
	3,   8,  12,  17,  25,  30,  41,  43,
	9,   11, 18,  24,  31,  40,  44,  53,
	10,  19, 23,  32,  39,  45,  52,  54,
	20,  22, 33,  38,  46,  51,  55,  60,
	21,  34, 37,  47,  50,  56,  59,  61,
	35,  36, 48,  49,  57,  58,  62,  63,
};

#if HUF_DEL_COUNT == 1
static unsigned long huffmanTotalDelay = 0;
static unsigned long huffmanDCDelay = 0;
static unsigned long huffmanACDelay = 0;
static unsigned long huffmanNonZeroCoef = 0;
static unsigned long huffmanAClengthsTotal = 0;
static unsigned long huffmanIsInHufCodesTotal = 0;
static unsigned long huffmanIsInHufCodesdelay = 0;
#endif
/***************************************************************************/


//struct stBlock
//{
//     int value;					// Decodes to.
//     int length;				// Length in bits.
//     unsigned short int code;	// 2 byte code (variable length)
//};

/***************************************************************************/


//struct stHuffmanTable
//{
//	unsigned char	m_length[17];		// 17 values from jpg file,
//										// k =1-16 ; L[k] indicates the number of Huffman codes of length k
//	unsigned char	m_hufVal[257];		// 256 codes read in from the jpeg file
//
//	int				m_numBlocks;
//	stBlock			m_blocks[1024];
//};


//struct stComponent
//{
//  unsigned int			m_hFactor;
//  unsigned int			m_vFactor;
//  float *				m_qTable;			// Pointer to the quantisation table to use
//  stHuffmanTable*		m_acTable;
//  stHuffmanTable*		m_dcTable;
//  short int				m_DCT[65];			// DCT coef
//  int					m_previousDC;
//};


//struct stJpegData
//{
//	unsigned char*		m_rgb;				// Final Red Green Blue pixel data
////	unsigned char		m_rgb[RGB_PIX_SIZE];		///
//	unsigned int		m_width;			// Width of image
//	unsigned int		m_height;			// Height of image
//
//	const unsigned char*m_stream;			// Pointer to the current stream
//	int					m_restart_interval;
//
//	stComponent			m_component_info[COMPONENTS];
//
//	float				m_Q_tables[COMPONENTS][64];	// quantization tables
//	stHuffmanTable		m_HTDC[HUFFMAN_TABLES];		// DC huffman tables
//	stHuffmanTable		m_HTAC[HUFFMAN_TABLES];		// AC huffman tables
//
//	// Temp space used after the IDCT to store each components
//	unsigned char		m_Y[64*4];
//	unsigned char		m_Cr[64];
//	unsigned char		m_Cb[64];
//
//	// Internal Pointer use for colorspace conversion, do not modify it !!!
////	unsigned char		m_colourspace[COLOURSPACE_SIZE];
//	unsigned char *		m_colourspace;
//};


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

/***************************************************************************/

// Clamp our integer between 0 and 255
inline unsigned char Clamp(int i)
{
	if (i<0)
		return 0;
	else if (i>255)
		return 255;
	else
		return i;
}

///***************************************************************************/
//
//#ifdef JPG_SPEEDUP
//
//#define GenHuffCodes( num_codes, arr, huffVal )				\
//{															\
//     int hufcounter = 0;									\
//     int codelengthcounter = 1;								\
//															\
//     for(int cc=0; cc< num_codes; cc++)						\
//     {														\
//		 while ( arr[cc].length > codelengthcounter )		\
//		 {													\
//			hufcounter = hufcounter << 1;					\
//			codelengthcounter++;							\
//		 }													\
//															\
//		 arr[cc].code = hufcounter;							\
//		 arr[cc].value = huffVal[cc];						\
//		 hufcounter = hufcounter + 1;						\
//	}														\
//}															\
//
//#else
//
//void GenHuffCodes( int num_codes, stBlock* arr, unsigned char* huffVal )
//{
//     int hufcounter = 0;
//     int codelengthcounter = 1;
//
//
//     for(int cc=0; cc< num_codes; cc++)
//     {
//		 while ( arr[cc].length > codelengthcounter )
//		 {
//			hufcounter = hufcounter << 1;
//			codelengthcounter++;
//		 }
//
//		 arr[cc].code = hufcounter;
//		 arr[cc].value = huffVal[cc];
//		 hufcounter = hufcounter + 1;
//	}
//}
//
//#endif //JPG_SPEEDUP


/***************************************************************************/

float C(int u)
{
    if (u == 0)
         return (1.0f/sqrtf(2));
    else
         return 1.0f;
}


//int func(int x, int y, const int block[8][8])
int IDCT_calc(int x, int y, short int block[8][8])
{
//#pragma HLS ARRAY_PARTITION variable=block complete dim=1
	const float PI = 3.14f;
    float sum=0;
//    float cos_coef_x = (2*x+1) * PI / 16;
//    float cos_coef_y = (2*y+1) * PI / 16;
    IDCT_calc_vertical_loop:
    for( int u=0; u<8; u++)
    {
//#pragma HLS UNROLL
    	IDCT_calc_horizontal_loop:
         for(int v=0; v<8; v++)
         {
//#pragma HLS UNROLL
//#pragma HLS PIPELINE II=1
             sum += ( C(u) * C(v) ) * block[u][v] * cosf( ((2*x+1) * u * PI) / 16)  * cosf( ((2*y+1) * v * PI) / 16);
//             sum += ( C(u) * C(v) ) * block[u][v] * cosf( cos_coef_x * u )  * cosf( cos_coef_y * v );
//             sum += C(v);
         }
    }
//    return (int) ((1.0/4.0) * sum);
    return (int) (0.25 * sum);
}

//void PerformIDCT(int outBlock[8][8], const int inBlock[8][8])
void PerformIDCT(short int outBlock[8][8], short int inBlock[8][8])
{
	PerformIDCT_vertical_loop:
	for(int y=0; y<8; y++)
	{
//#pragma HLS UNROLL
		PerformIDCT_horizontal_loop:
		for(int x=0; x<8; x++)
		{
//#pragma HLS UNROLL
			outBlock[x][y]  =  IDCT_calc( x, y, inBlock);
			outBlock[x][y] += 128;
		}
	}
}

/***************************************************************************/

//void DequantizeBlock( int block[64], const float quantBlock[64] )
void DequantizeBlock( short int block[64], short int coef[64], float quantBlock[64] )
{
	DequantizeBlock_loop:
	for( int c=0; c<64; c++)
	{
//#pragma HLS UNROLL
		block[c] = (int)(coef[c] * quantBlock[c]);
	}
}

/***************************************************************************/

//void DeZigZag(int outBlock[64], const int inBlock[64])
void DeZigZag(short int outBlock[64], short int inBlock[64])
{
	DeZigZag_loop:
	for(int i=0; i<64; i++)
	{
		outBlock[ i ] = inBlock[ZigZagArray[i]];
	}
}

/***************************************************************************/

//void TransformArray(int outArray[8][8], const int inArray[64])

void TransformArray(short int outArray[8][8], short  int inArray[64])
{
	int cc = 0;
	TransformArray_vertical_loop:
	for( int y=0; y<8; y++)
	{
		TransformArray_horizontal_loop:
		for( int x=0; x<8; x++)
		{
			outArray[x][y]  =  inArray[cc];
			cc++;
		}
	}
}

/***************************************************************************/

void DumpDecodedBlock(int val[8][8])
{
	dprintf("# Decoded 8x8 Block#\n");
	for( int y=0; y<8; y++)
	{
		for( int x=0; x<8; x++)
		{
			dprintf("%2x ", val[x][y]);
		}
		dprintf("\n");
	}
}

/***************************************************************************/

inline void DecodeSingleBlock(short int compDCT[64], float compqTable[64], unsigned char outputBuf[64*4], unsigned char stride)
//inline void DecodeSingleBlock(short int compDCT[65], float compqTable[64], unsigned char *outputBuf, int stride)
{
//	short* inptr    = compDCT;
//	float* quantptr = comp->m_qTable;


	// Create a temp 8x8, i.e. 64 array for the data
//	int data[64] = {0};

	short int coef[64];

	// Copy our data into the temp array
	Decode_SingleBlock_data_cp_loop:
	for (int i=0; i<64; i++)
	{

		coef[i] = compDCT[i];
//		data[i] = inptr[i];
//		printf("data[i]  %d\n",data[i]);
	}

	// De-Quantize
	short int data[64];
	DequantizeBlock(data, coef, compqTable);

//	for (int i=0; i<64; i++)
//	{
//		printf("new data[i]  %d\n",data[i]);
//	}

	// De-Zig-Zag
//	int block[64] = {0};
//	short int block[64] = {0};
	short int block[64];
	DeZigZag(block, data);

	// Create an 8x8 array
//	int arrayBlock[8][8]={0};
//	short int arrayBlock[8][8]={0};
	short int arrayBlock[8][8];
	TransformArray(arrayBlock, block);

	// Inverse DCT
//	int val[8][8]={0};
//	short int val[8][8]={0};
	short int val[8][8];
	PerformIDCT(val, arrayBlock);

	// Level Shift each element (i.e. add 128), and copy to our
	// outputf
//	unsigned char *outptr = outputBuf;
	SingleBlock_Save_vertical_loop:
	for (int y = 0; y < 8; y++)
	{
		SingleBlock_Save_horizontal_loop:
		for (int x=0; x<8; x++)
		{
//			val[x][y] += 128;

//			outptr[x] = Clamp(val[x][y]);
			outputBuf[x] = Clamp(val[x][y]);
		}

//		outptr += stride;
		outputBuf += stride;
	}

	#ifndef JPG_SPEEDUP
	DumpDecodedBlock(val);
	#endif
}

//inline void DecodeSingleBlock(stComponent *comp, unsigned char *outputBuf, int stride)
//{
//	short* inptr    = comp->m_DCT;
////	float* quantptr = comp->m_qTable;
//
//
//	// Create a temp 8x8, i.e. 64 array for the data
//	int data[64] = {0};
//
//	// Copy our data into the temp array
//	for (int i=0; i<64; i++)
//	{
//		data[i] = inptr[i];
////		printf("data[i]  %d\n",data[i]);
//	}
//
//	// De-Quantize
////	DequantizeBlock(data, quantptr);
//	DequantizeBlock(data, comp->m_qTable);
//
////	for (int i=0; i<64; i++)
////	{
////		printf("new data[i]  %d\n",data[i]);
////	}
//
//	// De-Zig-Zag
//	int block[64] = {0};
//	DeZigZag(block, data);
//
//	// Create an 8x8 array
//	int arrayBlock[8][8]={0};
//	TransformArray(arrayBlock, block);
//
//	// Inverse DCT
//	int val[8][8]={0};
//	PerformIDCT(val, arrayBlock);
//
//	// Level Shift each element (i.e. add 128), and copy to our
//	// output
//	unsigned char *outptr = outputBuf;
//	for (int y = 0; y < 8; y++)
//	{
//		for (int x=0; x<8; x++)
//		{
//			val[x][y] += 128;
//
//			outptr[x] = Clamp(val[x][y]);
//		}
//
//		outptr += stride;
//	}
//
//	#ifndef JPG_SPEEDUP
//	DumpDecodedBlock(val);
//	#endif
//}

/***************************************************************************/

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

///***************************************************************************/
//
//// Takes two array of bits, and build the huffman table for size, and code
//
///***************************************************************************/
//inline void BuildHuffmanTable(const unsigned char *bits, const unsigned char *stream, stHuffmanTable *HT)
//{
//	for (int j=1; j<=16; j++)
//	{
//		HT->m_length[j] = bits[j];
//	}
//
//	// Work out the total number of codes
//	int numBlocks = 0;
//	for (int i=1; i<=16; i++)
//	{
//		numBlocks += HT->m_length[i];
//	}
//	HT->m_numBlocks = numBlocks;
//
//	// Fill in the data our our blocks, so we know how many bits each
//	// one is
//	int c=0;
//	for (int i=1; i<=16; i++)
//	{
//		for (int j=0; j<HT->m_length[i]; j++)
//		{
//			HT->m_blocks[c].length = i;
//			c++;
//		}
//
//	}
//
//	GenHuffCodes(HT->m_numBlocks, HT->m_blocks, HT->m_hufVal);
//}
//
///***************************************************************************/
//
//inline void PrintSOF(const unsigned char *stream)
//{
//	int width;
//	int height;
//	int nr_components;
//	int precision;
//
//	const char *nr_components_to_string[] =	{	"????",
//												"Grayscale",
//												"????",
//												"YCbCr",
//												"CYMK" };
//
//	precision = stream[2];
//	height = BYTE_TO_WORD(stream+3);
//	width  = BYTE_TO_WORD(stream+5);
//	nr_components = stream[7];
//
//	dprintf("> SOF marker\n");
//	dprintf("Size:%dx%d nr_components:%d (%s)  precision:%d\n",
//							width, height,
//							nr_components,
//							nr_components_to_string[nr_components],
//							precision);
//}
//
///***************************************************************************/
//
////inline
//int ParseSOF(stJpegData *jdata, const unsigned char *stream)
//{
//	/*
//	SOF		16		0xffc0		Start Of Frame
//	Lf		16		3Nf+8		Frame header length
//	P		8		8			Sample precision
//	Y		16		0-65535		Number of lines
//	X		16		1-65535		Samples per line
//	Nf		8		1-255		Number of image components (e.g. Y, U and V).
//
//	---------Repeats for the number of components (e.g. Nf)-----------------
//	Ci		8		0-255		Component identifier
//	Hi		4		1-4			Horizontal Sampling Factor
//	Vi		4		1-4			Vertical Sampling Factor
//	Tqi		8		0-3			Quantization Table Selector.
//	*/
//
//	PrintSOF(stream);
//
//	int height = BYTE_TO_WORD(stream+3);
//	int width  = BYTE_TO_WORD(stream+5);
//	int nr_components = stream[7];
//
//	stream += 8;
//	for (int i=0; i<nr_components; i++)
//	{
//		int cid				= *stream++;
//		int sampling_factor = *stream++;
//		int Q_table			= *stream++;
//
//		stComponent *c = &jdata->m_component_info[cid];
//		c->m_vFactor = sampling_factor&0xf;
//		c->m_hFactor = sampling_factor>>4;
//		c->m_qTable = jdata->m_Q_tables[Q_table];
//
//		dprintf("Component:%d  factor:%dx%d  Quantization table:%d\n",
//				cid,/mnt/DATA/ANASTAS/Vivado_lin/jpeg_prj
//				c->m_vFactor,
//				c->m_hFactor,
//				Q_table );
//	}
//	jdata->m_width = width;
//	jdata->m_height = height;
//
//	return 0;
//}

///***************************************************************************/
//
//inline void BuildQuantizationTable(float *qtable, const unsigned char *ref_table)
//{
//	int c = 0;
//
//	for (int i=0; i<8; i++)
//	{
//		for (int j=0; j<8; j++)
//		{
//			unsigned char val = ref_table[c];
//
//			qtable[c] = val;
//			c++;
//		}
//	}
//}


///***************************************************************************/
//
//inline int ParseDQT(stJpegData *jdata, const unsigned char *stream)
//{
//	int length, qi;
//	float *table;
//
//	dprintf("> DQT marker\n");
//	length = BYTE_TO_WORD(stream) - 2;
//	stream += 2;	// Skip length
//
//	while (length>0)
//	{
//		qi = *stream++;
//
//		int qprecision = qi>>4;	 // upper 4 bits specify the precision
//		int qindex     = qi&0xf; // index is lower 4 bits
//
//		if (qprecision)
//		{
//			// precision in this case is either 0 or 1 and indicates the precision
//			// of the quantized values;
//			// 8-bit (baseline) for 0 and  up to 16-bit for 1
//
//			dprintf("Error - 16 bits quantization table is not supported\n");
//			///DBG_HALT;
//		}
//
//		if (qindex>=4)
//		{
//			dprintf("Error - No more 4 quantization table is supported (got %d)\n", qi);
//			///DBG_HALT;
//		}
//
//		// The quantization table is the next 64 bytes
//		table = jdata->m_Q_tables[qindex];
//
//		// the quantization tables are stored in zigzag format, so we
//		// use this functino to read them all in and de-zig zag them
//		BuildQuantizationTable(table, stream);
//		stream += 64;
//		length -= 65;
//	}
//	return 0;
//}
//
///***************************************************************************/
//
//inline int ParseSOS(stJpegData *jdata, const unsigned char *stream)
//{
//	/*
//	SOS		16		0xffd8			Start Of Scan
//	Ls		16		2Ns + 6			Scan header length
//	Ns		8		1-4				Number of image components
//	Csj		8		0-255			Scan Component Selector
//	Tdj		4		0-1				DC Coding Table Selector
//	Taj		4		0-1				AC Coding Table Selector
//	Ss		8		0				Start of spectral selection
//	Se		8		63				End of spectral selection
//	Ah		4		0				Successive Approximation Bit High
//	Ai		4		0				Successive Approximation Bit Low
//	*/
//
//	unsigned int nr_components = stream[2];
//
//	dprintf("> SOS marker\n");
//
//	if (nr_components != 3)
//	{
//		dprintf("Error - We only support YCbCr image\n");
//		///DBG_HALT;
//	}
//
//
//	stream += 3;
//	for (unsigned int i=0;i<nr_components;i++)
//	{
//		unsigned int cid   = *stream++;
//		unsigned int table = *stream++;
//
//		if ((table&0xf)>=4)
//		{
//			dprintf("Error - We do not support more than 2 AC Huffman table\n");
//			///DBG_HALT;
//		}
//		if ((table>>4)>=4)
//		{
//			dprintf("Error - We do not support more than 2 DC Huffman table\n");
//			///DBG_HALT;
//		}
//		dprintf("ComponentId:%d  tableAC:%d tableDC:%d\n", cid, table&0xf, table>>4);
//
//		jdata->m_component_info[cid].m_acTable = &jdata->m_HTAC[table&0xf];
//		jdata->m_component_info[cid].m_dcTable = &jdata->m_HTDC[table>>4];
//	}
//	jdata->m_stream = stream+3;
//	return 0;
//}
//
///***************************************************************************/
//
//inline int ParseDHT(stJpegData *jdata, const unsigned char *stream)
//{
//	/*
//	u8 0xff
//	u8 0xc4 (type of segment)
//	u16 be length of segment
//	4-bits class (0 is DC, 1 is AC, more on this later)
//	4-bits table id
//	array of 16 u8 number of elements for each of 16 depths
//	array of u8 elements, in order of depth
//	*/
//
//	unsigned int count, i;
//	unsigned char huff_bits[17];
//	int length, index;
//
//	length = BYTE_TO_WORD(stream) - 2;
//	stream += 2;	// Skip length
//
//	dprintf("> DHT marker (length=%d)\n", length);
//
//	while (length>0)
//	{
//		index = *stream++;
//
//		// We need to calculate the number of bytes 'vals' will takes
//		huff_bits[0] = 0;
//		count = 0;
//		for (i=1; i<17; i++)
//		{
//			huff_bits[i] = *stream++;
//			count += huff_bits[i];
//		}
//
//		if (count > 256)
//		{
//			dprintf("Error - No more than 1024 bytes is allowed to describe a huffman table");
//			///DBG_HALT;
//		}
//		if ( (index &0xf) >= HUFFMAN_TABLES)
//		{
//			dprintf("Error - No mode than %d Huffman tables is supported\n", HUFFMAN_TABLES);
//			///DBG_HALT;
//		}
//		dprintf("Huffman table %s n%d\n", (index&0xf0)?"AC":"DC", index&0xf);
//		dprintf("Length of the table: %d\n", count);
//
//		if (index & 0xf0 )
//		{
//			unsigned char* huffval = jdata->m_HTAC[index&0xf].m_hufVal;
//			for (i = 0; i < count; i++)
//				huffval[i] = *stream++;
//
//			BuildHuffmanTable(huff_bits, stream, &jdata->m_HTAC[index&0xf]); // AC
//		}
//		else
//		{
//			unsigned char* huffval = jdata->m_HTDC[index&0xf].m_hufVal;
//			for (i = 0; i < count; i++)
//				huffval[i] = *stream++;
//
//			BuildHuffmanTable(huff_bits, stream, &jdata->m_HTDC[index&0xf]); // DC
//		}
//
//		length -= 1;
//		length -= 16;
//		length -= count;
//	}
//	dprintf("< DHT marker\n");
//	return 0;
//}
//
///***************************************************************************/
//
//inline int ParseJFIF(stJpegData *jdata, const unsigned char *stream)
//{
//	int chuck_len;
//	int marker;
//	int sos_marker_found = 0;
//	int dht_marker_found = 0;
//
//	// Parse marker
//	while (!sos_marker_found)
//	{
//		if (*stream++ != 0xff)
//		{
//			goto bogus_jpeg_format;
//		}
//
//		// Skip any padding ff byte (this is normal)
//		while (*stream == 0xff)
//		{
//			stream++;
//		}
//
//		marker = *stream++;
//		chuck_len = BYTE_TO_WORD(stream);
//
//		switch (marker)
//		{
//			case SOF:
//			{
//				if (ParseSOF(jdata, stream) < 0)
//					return -1;
//			}
//			break;
//
//			case DQT:
//			{
//				if (ParseDQT(jdata, stream) < 0)
//					return -1;
//			}
//			break;
//
//			case SOS:
//			{
//				if (ParseSOS(jdata, stream) < 0)
//					return -1;
//				sos_marker_found = 1;
//			}
//			break;
//
//			case DHT:
//			{
//				if (ParseDHT(jdata, stream) < 0)
//					return -1;
//				dht_marker_found = 1;
//			}
//			break;
//
//			// The reason I added these additional skips here, is because for
//			// certain jpg compressions, like swf, it splits the encoding
//			// and image data with SOI & EOI extra tags, so we need to skip
//			// over them here and decode the whole image
//			case SOI:
//			case EOI:
//			{
//				chuck_len = 0;
//				break;
//			}
//			break;
//
//			case 0xDD: //DRI: Restart_markers=1;
//			{
//				jdata->m_restart_interval = BYTE_TO_WORD(stream);
//				dprintf("DRI - Restart_marker\n");
//			}
//			break;
//
//			case APP0:
//			{
//				dprintf("APP0 Chunk ('txt' information) skipping\n");
//			}
//			break;
//
//			default:
//			{
//				dprintf("ERROR> Unknown marker %2.2x\n", marker);
//			}
//			break;
//		}
//
//		stream += chuck_len;
//	}
//
//	if (!dht_marker_found)
//	{
//		dprintf("ERROR> No Huffman table loaded\n");
//		///DBG_HALT;
//	}
//
//	return 0;
//
//	bogus_jpeg_format:
//	dprintf("ERROR> Bogus jpeg format\n");
//	///DBG_HALT;
//	return -1;
//}

///***************************************************************************/
//
////inline int JpegParseHeader(stJpegData *jdata, const unsigned char *buf, unsigned int size)
//inline int JpegParseHeader(stJpegData *jdata, const unsigned char *buf, unsigned int size)
//{
//	// Identify the file
//	if ((buf[0] != 0xFF) || (buf[1] != SOI))
//	{
//		dprintf("Not a JPG file ?\n");
//		///DBG_HALT;
//		return -1;
//	}
//
//	const unsigned char* startStream = buf+2;
//	const int fileSize = size-2;
//
//	dprintf("-|- File thinks its size is: %d bytes\n", fileSize);
//
//	int ret = ParseJFIF(jdata, startStream);
//
//
////	printf("\n ret: %d\n", ret);
//
//	return ret;
//}

/***************************************************************************/

/////inline
//void JpegGetImageSize(stJpegData *jdata, unsigned int *width, unsigned int *height)
//{
//	*width  = jdata->m_width;
//	*height = jdata->m_height;
//}

/***************************************************************************/

//unsigned int g_reservoir = 0;
//unsigned int g_nbits_in_reservoir = 0;

#ifdef JPG_SPEEDUP
//#define FillNBits(stream, nbits_wanted)						\
//{															\
//	while ((int)g_nbits_in_reservoir<nbits_wanted)			\
//	{														\
//		unsigned char c = *(stream)++;						\
//		printf("c size:%d\n\n",c);						\
//		printf("stream size:%d\n\n",*(stream));						\
//		g_reservoir <<= 8;									\
//		if (c == 0xff && (*stream) == 0x00)				\
//		{													\
//		/*	(stream)++;		*/							\
//		}													\
//		g_reservoir |= c;									\
//		g_nbits_in_reservoir+=8;							\
//	}														\
//}
#if (STREAM_POINTER != 0)
//#define FillNBits(stream, nbits_wanted)						\
//{															\
//	int lll;			\
//	while ((int)g_nbits_in_reservoir<nbits_wanted)			\
//	{														\
//		/**printf("*(*stream) is:%d\n",*(*stream));	*/		\
//		unsigned char c = *(*stream);						\
//		(*stream)+=1;						\
//		ll++;									\
//		/*c++;		*/				\
//		/*printf("*(*stream)++ is:%d\n",c);		*/				\
//		/*printf("*(*stream) is:%d\n",*(*stream));		*/				\
//		/*printf("(**stream) is:%d\n",(**stream));		*/				\
//printf("(*stream):%p\n",(stream));							\
///*printf("stream is:%d\n",stream);		*/			\
//		g_reservoir <<= 8;									\
//		if (c == 0xff && (**stream) == 0x00)				\
//		{													\
//			(*stream)+=1;									\
//			ll++;									\
//		}													\
//		/*printf("*(*stream) is:%d\n",*(*stream));	*/		\
//		printf("ll is:%d\n",ll);						\
//		g_reservoir |= c;									\
//		g_nbits_in_reservoir+=8;							\
//	}														\
//	/*printf("end of while\n");			*/		\
//}
//inline
void  FillNBits(unsigned char **stream, int nbits_wanted)						\
{
	int lll;			\
	while ((int)g_nbits_in_reservoir<nbits_wanted)
	{
		/**printf("*(*stream) is:%d\n",*(*stream));	*/
		unsigned char c = *(*stream);
		(*stream)+=1;						\
		ll++;									\
		/*c++;		*/				\
		/*printf("*(*stream)++ is:%d\n",c);		*/
		/*printf("*(*stream) is:%d\n",*(*stream));		*/				\
		/*printf("(**stream) is:%d\n",(**stream));		*/				\
//printf("(*stream):%p\n",(*stream));
/*printf("stream is:%d\n",stream);		*/			\
		g_reservoir <<= 8;									\
		if (c == 0xff && (**stream) == 0x00)				\
		{													\
			(*stream)+=1;									\
			ll++;									\
		}													\
		/*printf("*(*stream) is:%d\n",*(*stream));	*/		\
		//printf("ll is:%d\n",ll);
		g_reservoir |= c;									\
		g_nbits_in_reservoir+=8;							\
	}														\
	/*printf("end of while\n");			*/		\
}
#else
//inline
//void  FillNBits(unsigned char **stream, int nbits_wanted)						\
//{
//	int lll;			\
//	while ((int)g_nbits_in_reservoir<nbits_wanted)
//	{
//		/**printf("*(*stream) is:%d\n",*(*stream));	*/
//		unsigned char c = *(*stream);
//		(*stream)+=1;						\
//		ll++;									\
//		/*c++;		*/				\
//		/*printf("*(*stream)++ is:%d\n",c);		*/
//		/*printf("*(*stream) is:%d\n",*(*stream));		*/				\
//		/*printf("(**stream) is:%d\n",(**stream));		*/				\
////printf("(*stream):%p\n",(*stream));
///*printf("stream is:%d\n",stream);		*/			\
//		g_reservoir <<= 8;									\
//		if (c == 0xff && (**stream) == 0x00)				\
//		{													\
//			(*stream)+=1;									\
//			ll++;									\
//		}													\
//		/*printf("*(*stream) is:%d\n",*(*stream));	*/		\
//		//printf("ll is:%d\n",ll);
//		g_reservoir |= c;									\
//		g_nbits_in_reservoir+=8;							\
//	}														\
//	/*printf("end of while\n");			*/		\
//}
//#define FillNBits(stream, nbits_wanted)						\
//{															\
//	int lll;			\
//	while ((int)g_nbits_in_reservoir<nbits_wanted)			\
//	{														\
//		/**printf("*(*stream) is:%d\n",*(*stream));	*/		\
//		unsigned char c = *(*stream);						\
//		(*stream)+=1;						\
//		ll++;									\
//		/*c++;		*/				\
//		/*printf("*(*stream)++ is:%d\n",c);		*/				\
//		/*printf("*(*stream) is:%d\n",*(*stream));		*/				\
//		/*printf("(**stream) is:%d\n",(**stream));		*/				\
//printf("(*stream):%p\n",(stream));							\
///*printf("stream is:%d\n",stream);		*/			\
//		g_reservoir <<= 8;									\
//		if (c == 0xff && (**stream) == 0x00)				\
//		{													\
//			(*stream)+=1;									\
//			ll++;									\
//		}													\
//		/*printf("*(*stream) is:%d\n",*(*stream));	*/		\
//		printf("ll is:%d\n",ll);						\
//		g_reservoir |= c;									\
//		g_nbits_in_reservoir+=8;							\
//	}														\
//	/*printf("end of while\n");			*/		\
//}
//		#define FillNBits(stream, nbits_wanted)						\
//		{															\
//			int lll; \
//			for(;((int)g_nbits_in_reservoir<nbits_wanted); )\
//			{														\
//				/*lll++;		*/				\
//				unsigned char c = stream[0];						\
//				stream+=1;												\
//				g_reservoir <<= 8;									\
//				if (c == 0xff && stream[0] == 0x00)				\
//				{													\
//				/*	(*stream)++;	*/ stream+=1;								\
//				}													\
//				printf("stream:%d\n",c);						\
//				printf("stream+1 is:%d\n",stream[0]);						\
//				printf("stream+2 size:%d\n",stream[0+1]);						\
//				printf("l is:%d\t\t\t\t\n",lll);						\
//				g_reservoir |= c;									\
//				g_nbits_in_reservoir+=8;							\
//			}														\
//				/*printf("end of while\n");	*				\
//		}
//#define FillNBits(stream, nbits_wanted)						\
//{															\
//	int lll;			\
//	while ((int)g_nbits_in_reservoir<nbits_wanted)			\
//	{														\
//		printf("*(*stream) is:%d\n",(*stream));			\
//		unsigned char c = *stream;						\
//		stream+=2;						\
//		ll++;									\
//		/*c++;		*/				\
//		/*printf("*(*stream)++ is:%d\n",c);		*/				\
//		printf("*(*stream) is:%d\n",*stream);						\
//		/*printf("(**stream) is:%d\n",(**stream));		*/				\
///*printf("(*stream):%d\n",(*stream));			*/				\
///*printf("stream is:%d\n",stream);		*/			\
//		g_reservoir <<= 8;									\
//		if (c == 0xff && *stream == 0x00)				\
//		{													\
//			stream+=1;									\
//			ll++;									\
//		}													\
//		/*printf("*(*stream) is:%d\n",*(*stream));	*/		\
//		printf("ll is:%d\n",ll);						\
//		g_reservoir |= c;									\
//		g_nbits_in_reservoir+=8;							\
//	}														\
//	/*printf("end of while\n");		*/			\
//}
#endif

//#if (STREAM_POINTER != 0)
////inline short GetNBits(const unsigned char** stream, int nbits_wanted)
//inline short GetNBits(unsigned char** stream, int nbits_wanted)
//#else
//inline short GetNBits(unsigned char* stream, int nbits_wanted)
////inline short GetNBits(const unsigned char* stream, int nbits_wanted)
//#endif
//{
//	/********************** FillNBits ***************************/\
//		while ((int)g_nbits_in_reservoir<k)			\
//		{														\
//			unsigned char c = jdata->m_stream[jdata->stream_index + 0];						\
//			jdata->stream_index++;						\
//			printf("c size:%d\n\n",c);						\
//			printf("stream size:%d\n\n",jdata->m_stream[jdata->stream_index + 0]);						\
//			g_reservoir <<= 8;									\
//			if (c == 0xff && jdata->m_stream[jdata->stream_index + 0] == 0x00)				\
//			{													\
//				jdata->stream_index++;						\
//			}													\
//			g_reservoir |= c;									\
//			g_nbits_in_reservoir+=8;							\
//		}	\
////	FillNBits(stream, nbits_wanted);
//
//	short result = ((g_reservoir)>>(g_nbits_in_reservoir-(nbits_wanted)));
//
//	g_nbits_in_reservoir -= (nbits_wanted);
//	g_reservoir &= ((1U<<g_nbits_in_reservoir)-1);
//
//	return result;
//}

//#if (STREAM_POINTER != 0)
////inline int LookNBits(const unsigned char** stream, int nbits_wanted)
//inline int LookNBits(unsigned char** stream, int nbits_wanted)
//#else
////inline int LookNBits(const unsigned char* stream, int nbits_wanted)
//inline int LookNBits(unsigned char* stream, int nbits_wanted)
//#endif
//{
//	/********************** FillNBits ***************************/
//		while ((int)g_nbits_in_reservoir<nbits_wanted)			\
//		{														\
//			unsigned char c = *(stream)++;						\
//			printf("c size:%d\n\n",c);						\
//			printf("stream size:%d\n\n",*(stream));						\
//			g_reservoir <<= 8;									\
//			if (c == 0xff && (*stream) == 0x00)				\
//			{													\
//			/*	(stream)++;		*/							\
//			}													\
//			g_reservoir |= c;									\
//			g_nbits_in_reservoir+=8;							\
//		}
////	FillNBits(stream, nbits_wanted);
//
//	int result = ((g_reservoir)>>(g_nbits_in_reservoir-(nbits_wanted)));
//	return result;
//}

//#define SkipNBits(stream, nbits_wanted)						\
//{															\
//	/********************** FillNBits ***************************/\
//		while ((int)g_nbits_in_reservoir<nbits_wanted)			\
//		{														\
//			unsigned char c = *(stream);						\
//			printf("c size:%d\n\n",c);						\
//			printf("stream size:%d\n\n",*(stream));						\
//			g_reservoir <<= 8;									\
//			if (c == 0xff && (*stream) == 0x00)				\
//			{													\
//			/*	(stream)++;		*/							\
//			}													\
//			g_reservoir |= c;									\
//			g_nbits_in_reservoir+=8;							\
//		}	\
///*	FillNBits(stream, nbits_wanted);	*/					\
//															\
//	g_nbits_in_reservoir -= (nbits_wanted);					\
//	g_reservoir &= ((1U<<g_nbits_in_reservoir)-1);			\
//}															\


#else

inline void FillNBits(const unsigned char** stream, int& nbits_wanted)
{
	while ((int)g_nbits_in_reservoir<nbits_wanted)
	{
		/*const */unsigned char c = *(*stream)++;
		g_reservoir <<= 8;
		//if (c == 0xff && (**stream) != 0x00)
		//{
		//	//DBG_HALT;
		//	*(*stream)++;
		//	*(*stream)++;
		//	c = *(*stream)++;

		//	g_reservoir = 0;
		//	g_nbits_in_reservoir = 0;
		//}

		if (c == 0xff && (**stream) == 0x00)
		{
			(*stream)++;
		}
		g_reservoir |= c;
		g_nbits_in_reservoir+=8;
	}
}

inline short GetNBits(const unsigned char** stream, int nbits_wanted)
{
	FillNBits(stream, nbits_wanted);

	short result = ((g_reservoir)>>(g_nbits_in_reservoir-(nbits_wanted)));

	g_nbits_in_reservoir -= (nbits_wanted);
	g_reservoir &= ((1U<<g_nbits_in_reservoir)-1);

	/*
	// Could do the sign conversion here!
	if (result < (short)(1UL<<((nbits_wanted)-1)))
	{
		result = result + (short)(0xFFFFFFFFUL<<(nbits_wanted))+1;
	}
	*/
	return result;
}

inline int LookNBits(const unsigned char** stream, int nbits_wanted)
{
	FillNBits(stream, nbits_wanted);

	int result = ((g_reservoir)>>(g_nbits_in_reservoir-(nbits_wanted)));
	return result;
}

inline void SkipNBits(const unsigned char** stream, int& nbits_wanted)
{
//	FillNBits(stream, nbits_wanted);
//
//	g_nbits_in_reservoir -= (nbits_wanted);
//	g_reservoir &= ((1U<<g_nbits_in_reservoir)-1);
}

#endif //JPG_SPEEDUP

/***************************************************************************/


bool IsInHuffmanCodes(int code, int numCodeBits, int numBlocks, stBlock *blocks, int* outValue)
{
#if SET_CONST_VALUES != 0
	numBlocks = 162;
#endif

//	printf("numBlocks: %d\n\n",numBlocks);

//	 hls::stream<stBlock> *blocks,
//	#pragma HLS INTERFACE port=blocks axis
//	#pragma HLS data_pack variable=blocks	// force all elements of struct to be packed in together
	// dindt work??
	bool returnResult = false;

//	for (int j=0; j<162; j++)
	IsinHuffmanCodes_loop:
#if HUF_DEL_COUNT == 1
	uint8_t j;
	IsinHuffmanCodes_loop:
	for (j=0; j<numBlocks; j++)
#else
	for (uint8_t j=0; j<numBlocks; j++)
#endif
	{
#if PIC == 1
#pragma HLS LOOP_TRIPCOUNT min=5 max=5 avg=5
#elif PIC == 2
#pragma HLS LOOP_TRIPCOUNT min=2 max=2 avg=2
#elif PIC == 3
#pragma HLS LOOP_TRIPCOUNT min=4 max=4 avg=4
#else
#pragma HLS LOOP_TRIPCOUNT min=1 max=162 avg=65
#endif
		int hufhCode		= blocks[j].code;
		int hufCodeLenBits	= blocks[j].length;
		int hufValue		= blocks[j].value;

		// We've got a match!
		if ((code==hufhCode) && (numCodeBits==hufCodeLenBits))
		{
			*outValue = hufValue;
//			return true;
			returnResult = true;
//			j = numBlocks;	// break loop

#if HUF_DEBUG == 1
			printf("%d\t. index of Huf Table\n",j);
#endif
		}
	}
#if HUF_DEL_COUNT == 1
			huffmanIsInHufCodesdelay += j+1;
			huffmanIsInHufCodesTotal += j+1;
#endif
//	return false;
	return returnResult;
}

/***************************************************************************/

int DetermineSign(int val, int nBits)
{
	bool negative = val < (1<<(nBits-1));

	if (negative)
	{
		// (-1 << (s)), makes the last bit a 1, so we have 1000,0000 for example for 8 bits

		val = val + (-1 << (nBits)) + 1;
	}

	// Else its unsigned, just return
	return val;
}

/***************************************************************************/

char g_bigBuf[1024] = {0};
char* IntToBinary(int val, int bits)
{
	for (int i=0; i<32; i++) g_bigBuf[i]='\0';

	int c = 0;
	for (int i=bits-1; i>=0; i--)
	{
		bool on = (val & (1<<i)) ? 1 : 0;
		g_bigBuf[c] = on ? '1' : '0';
		c++;
	}

	return &g_bigBuf[0];
}

/***************************************************************************/

void DumpHufCodes(stHuffmanTable* table)
{
	printf("HufCodes\n");
	printf("Num: %d\n", table->m_numBlocks);
	for (int i = 0; i<table->m_numBlocks; i++)
	{
		printf("%03d\t [%s]\n", i, IntToBinary(table->m_blocks[i].code, table->m_blocks[i].length));
	}
	printf("\n");

}

/***************************************************************************/
//
//void DumpDCTValues(short dct[64])
//{
//	dprintf("\n#Extracted DCT values from SOS#\n");
//	int c = 0;
//	for (int i=0; i<64; i++)
//	{
//		dprintf("% 4d  ", dct[c++]);
//
//		if ( (c>0) && (c%8==0) ) dprintf("\n");
//	}
//	dprintf("\n");
//}
//
//
//void ProcessRestart()
//{
//
//}

/***************************************************************************/
//void FillNBits(stJpegData *jdata, char limit)
void FillNBits(stHuffmanData *hufData, char limit)
{
	FillNBits_loop:
#if FILLNBITS_SPEEDUP == 1
	while(hufData->g_nbits_in_reservoir<limit)
#pragma HLS LOOP_TRIPCOUNT min=1 max=2
	{
//		unsigned char c = jdata->m_stream[jdata->stream_index];
//		jdata->stream_index++;
//		jdata->g_reservoir <<= 8;
//		jdata->g_reservoir |= c;
		hufData->g_reservoir = (hufData->g_reservoir << 8) | hufData->m_stream[hufData->stream_index];
		hufData->g_nbits_in_reservoir+=8;
		if (hufData->m_stream[hufData->stream_index] == 0xff && hufData->m_stream[hufData->stream_index+1] == 0x00)
		{
//			jdata->stream_index++;
			hufData->stream_index+=2;
		}
		else
			hufData->stream_index++;
	}
#else
	for(char byte_cnt=0; byte_cnt<2; byte_cnt++)
	{
		if(hufData->g_nbits_in_reservoir<limit)
		{
			unsigned char c = hufData->m_stream[hufData->stream_index];
			hufData->stream_index++;
			hufData->g_reservoir <<= 8;
			if (c == 0xff && hufData->m_stream[hufData->stream_index] == 0x00)
			{
				hufData->stream_index++;
			}
			hufData->g_reservoir |= c;
			hufData->g_nbits_in_reservoir+=8;
		}
	}
#endif
}

//void ProcessHuffmanBlock(stComponent *c, *stHuffmanTable HTDC, *stHuffmanTable HTAC,  int indx)
//void ProcessHuffmanBlock(stComponent *c, stHuffmanData *hufUnit,  int indx)
#if HUF_LUT_SPEEDUP == 1
#if HUF_DCT_SPEEDUP == 0
void ProcessHuffmanBlock(short int m_DCT[64], short int *previousDC, stHufLUT *m_hufLUT, stHuffmanData *hufUnit)
#else
void ProcessHuffmanBlock(short int DCT[64], short int *previousDC, stHufLUT *DC_hufLUT, stHufLUT *AC_hufLUT, stHuffmanData *hufUnit)
//void ProcessHuffmanBlock(short int DCT[64], short int *previousDC, stHufLUT DC_hufLUT[65536], stHufLUT AC_hufLUT[65536], stHuffmanData *hufUnit)
#endif
#else
#if HUF_DCT_SPEEDUP == 0
void ProcessHuffmanBlock(short int m_DCT[64], short int *previousDC, stHuffmanTable *HTDC, stHuffmanTable *HTAC, stHuffmanData *hufUnit)
#else
void ProcessHuffmanBlock(short int DCT[64], short int *previousDC, stHuffmanTable *HTDC, stHuffmanTable *HTAC, stHuffmanData *hufUnit)
#endif
#endif
//void ProcessHuffmanBlock(short int m_DCT[64], stHuffmanTable *m_HTDC, stHuffmanTable *m_HTAC, unsigned int *g_nbits_in_reservoir, unsigned int *g_reservoir, int *restart_interval)
//void ProcessHuffmanBlock(stJpegData* jdata,  int indx)
{
//	printf("gNbits->g_nbits_in_reservoir:%d\n\n",jdata->g_nbits_in_reservoir);

#if HUF_DEL_COUNT == 1
	huffmanDCDelay = huffmanACDelay = 0;
	huffmanIsInHufCodesdelay = 0;
#endif
	///DBG_ASSERT(indx>=0 && indx<4);
//	stComponent *c = &jdata->m_component_info[indx];

	// Start Huffman decoding

	// We memset it here, as later on we can just skip along, when we have lots
	// of leading zeros, for our AC run length encoding :)
#if HUF_DCT_SPEEDUP == 0
	short DCT_tcoeff[64];
	memset(DCT_tcoeff, 0, sizeof(DCT_tcoeff)); //Initialize DCT_tcoeff
#else
//	memset(DCT, 0, sizeof(DCT)); //Initialize DCT_tcoeff
	init_DCT_loop:
	for (char cnt=0; cnt<64; cnt++)
	{
		DCT[cnt] = 0;
	}
#endif

	bool found = false;
	int decodedValue = 0;

//
//////	DumpHufCodes(c->dcTable);
//////	DumpHufCodes(c->acTable);
////	DumpHufCodes(&jdata->m_Huffman.m_HTDC[indx]);
////	DumpHufCodes(&jdata->m_Huffman.m_HTAC[indx]);
//	DumpHufCodes(&hufUnit->m_HTDC[indx]);
//	DumpHufCodes(&hufUnit->m_HTAC[indx]);

//	printf("\nHuff Block:\n\n");

//	dprintf("\nHuff Block:\n\n");


	// Scan Decode Resync
	if (hufUnit->m_restart_interval>0)
	if ( hufUnit->m_stream[hufUnit->stream_index + 0]==0xff && hufUnit->m_stream[hufUnit->stream_index + 1]!=0x00)
	{
		// Something might be wrong, we should have had an interval marker set
		//DBG_ASSERT(jdata->m_restart_interval>0);
		printf("jdata->stream_index is:%d\n",hufUnit->stream_index);
//		jdata->m_stream += 2;
		hufUnit->stream_index += 2;
		hufUnit->g_reservoir = 0;
		hufUnit->g_nbits_in_reservoir = 0;

		// The value in the interval marker determines what number it will count
		// upto before looping back...i.e, for an interval marker of 4, we will
		// have 0xFFD0, 0xFFD1, 0xFFD2, 0xFFD3, 0xFFD4, 0xFFD0, 0xFFD1..etc..looping
		// so each time we get a resync, it will count from 0 to 4 then restart :)
	}


//	printf("jdata->stream_index is:%d\n",jdata->stream_index);

//	ProcessHuffmanDC(jdata);

	// First thing is get the 1 DC coefficient at the start of our 64 element
	// block
	Huffman_getCoef_DC_loop:
#if HUF_DEL_COUNT == 1
	int k;
	for (k=1; k<16; k++)
#else
	for (int k=1; k<16; k++)
#endif
	{
#if PIC == 1
#pragma HLS LOOP_TRIPCOUNT min=5 max=5 avg=5
#elif PIC == 2
#pragma HLS LOOP_TRIPCOUNT min=2 max=2 avg=2
#elif PIC == 3
#pragma HLS LOOP_TRIPCOUNT min=4 max=4 avg=4
#else

#endif
		// Keep grabbing one bit at a time till we find one thats a huffman code
#if (STREAM_POINTER != 0)
		int code = LookNBits(&jdata->m_stream, k);
#else
		/********************** FillNBits ***************************/
			DC_LookNBits_FillNBits_loop:
			FillNBits(hufUnit,k);

			int code = ((hufUnit->g_reservoir)>>(hufUnit->g_nbits_in_reservoir-(k)));

#endif
		// Check if its one of our huffman codes
//		if (IsInHuffmanCodes(code, k,  c->m_dcTable->m_numBlocks, c->m_dcTable->m_blocks, &decodedValue))

//		if (IsInHuffmanCodes(code, k,  hufUnit->m_HTDC[c->dcTable_index].m_numBlocks, \
//				hufUnit->m_HTDC[c->dcTable_index].m_blocks, &decodedValue))
#if HUF_LUT_SPEEDUP == 1
		if ( DC_hufLUT[code].DClength == k)
		{
				decodedValue = DC_hufLUT[code].DCvalue;
#else
		if (IsInHuffmanCodes(code, k,  HTDC->m_numBlocks, HTDC->m_blocks, &decodedValue))
		{
#endif
			// Skip over the rest of the bits now.
#if (STREAM_POINTER != 0)
			SkipNBits(&jdata->m_stream, k);
#else
			/********************** FillNBits ***************************/\
				DC_SkipNBits_FillNBits_loop:

			FillNBits(hufUnit,k);
			hufUnit->g_nbits_in_reservoir -= (k);
			hufUnit->g_reservoir &= ((1U<<hufUnit->g_nbits_in_reservoir)-1);
//			SkipNBits(jdata->m_stream, k);
#endif

			found = true;

			// The decoded value is the number of bits we have to read in next
			int numDataBits = decodedValue;

			// We know the next k bits are for the actual data
			if (numDataBits==0)
			{
//				DCT_tcoeff[0] = c->m_previousDC;
#if HUF_DCT_SPEEDUP == 0
				DCT_tcoeff[0] = *previousDC;
#else
				DCT[0] = *previousDC;
#endif
			}
			else
			{

				if (hufUnit->m_restart_interval>0)
//					if ( jdata->m_stream[0]==0xff && jdata->m_stream[1]!=0x00)
						if ( hufUnit->m_stream[hufUnit->stream_index + 0]==0xff && hufUnit->m_stream[hufUnit->stream_index + 1]!=0x00)
				{
//					jdata->m_stream += 2;
							hufUnit->stream_index += 2;

							hufUnit->g_reservoir = 0;
							hufUnit->g_nbits_in_reservoir = 0;
				}
#if (STREAM_POINTER != 0)
				short data = GetNBits(&jdata->m_stream, numDataBits);
#else
//				printf("numDataBits%d\n\n",numDataBits);
				/********************** FillNBits ***************************/\
					DC_GetNBits_FillNBits_loop:
					FillNBits(hufUnit,numDataBits);

				short data = ((hufUnit->g_reservoir)>>(hufUnit->g_nbits_in_reservoir-(numDataBits)));

				hufUnit->g_nbits_in_reservoir -= (numDataBits);
				hufUnit->g_reservoir &= ((1U<<hufUnit->g_nbits_in_reservoir)-1);
//				short data = GetNBits(jdata->m_stream, numDataBits);
#endif

				data = DetermineSign(data, numDataBits);

////				DCT_tcoeff[0] = data + c->m_previousDC;
////				c->m_previousDC = DCT_tcoeff[0];
#if HUF_DCT_SPEEDUP == 0
				DCT_tcoeff[0] = data + *previousDC;
				*previousDC = DCT_tcoeff[0];
#else
				DCT[0] = data + *previousDC;
				*previousDC = DCT[0];
#endif
			}
#if HUF_DEBUG == 1
			printf("%d\tbits wide DC Huf Code\n",k);
			printf("%d\tvalue of DC Huf Table\n",numDataBits);
#if HUF_LUT_SPEEDUP == 0
			printf("%d\t m_numBlocks\n",HTDC->m_numBlocks);
#endif
#endif
			// Found so we can exit out
			break;
		}
	}

	if (!found)
	{
		printf("-|- ##ERROR## We have a *serious* error, unable to find huffman code\n");
//		dprintf("-|- ##ERROR## We have a *serious* error, unable to find huffman code\n");
		///DBG_HALT;
	}
#if HUF_DEL_COUNT == 1
//	huffmanDCDelay += (k-2) * 162 * ISINHUF_SEARCH_FACTOR + \
//						(huffmanIsInHufCodesdelay * ISINHUF_SEARCH_FACTOR)  + HUF_DC_SEARCH_LATENCY;
	huffmanDCDelay += (huffmanIsInHufCodesdelay * ISINHUF_SEARCH_FACTOR + ISINHUF_INIT)  + HUF_DC_SEARCH_LATENCY;
	huffmanIsInHufCodesdelay = 0;
#endif
	//if ( jdata->m_stream[0]==0xff && jdata->m_stream[1]!=0x00)
	//{
	//	jdata->m_stream += 2;
	//	g_reservoir = 0;
	//	g_nbits_in_reservoir = 0;
	//}


	// Second, the 63 AC coefficient
	int nr=1;
	bool EOB_found=false;
//	Huffman_getEOB_AC_loop:
//	for(int ac_coef_cnt=0; ac_coef_cnt<64; ac_coef_cnt++)
//	{
//		printf("nr%d\n\n",nr);
//	if ( (nr<=63) && (!EOB_found) )
//		{

	Huffman_getEOB_AC_loop:
	while ( (nr<=63)&&(!EOB_found) )
	{
#if PIC == 1
#pragma HLS LOOP_TRIPCOUNT min=40 max=40
#elif PIC == 2
#pragma HLS LOOP_TRIPCOUNT min=6 max=6
#elif PIC == 3
#pragma HLS LOOP_TRIPCOUNT min=35 max=35
#else
#pragma HLS LOOP_TRIPCOUNT min=1 max=63
#endif
		int k = 0;
		Huffman_getCoef_AC_loop:
		for (k=1; k<=16; k++)
		{
#if PIC == 1
#pragma HLS LOOP_TRIPCOUNT min=5 max=5 avg=5
#elif PIC == 2
#pragma HLS LOOP_TRIPCOUNT min=2 max=2 avg=2
#elif PIC == 3
#pragma HLS LOOP_TRIPCOUNT min=4 max=4 avg=4
#else

#endif
			// Keep grabbing one bit at a time till we find one thats a huffman code
#if (STREAM_POINTER != 0)
			int code = LookNBits(&jdata->m_stream, k);
#else
			/********************** FillNBits ***************************/\
			AC_LookNBits_FillNBits_loop:
			FillNBits(hufUnit,k);
//
			int code = ((hufUnit->g_reservoir)>>(hufUnit->g_nbits_in_reservoir-(k)));


//			int code = LookNBits(jdata->m_stream, k);
#endif


			// Check if its one of our huffman codes

//					if (IsInHuffmanCodes(code, k,  hufUnit->m_HTAC[c->acTable_index].m_numBlocks, \
//												   hufUnit->m_HTAC[c->acTable_index].m_blocks, &decodedValue))
#if HUF_LUT_SPEEDUP == 1
			if (AC_hufLUT[code].AClength == k)
			{
				decodedValue = AC_hufLUT[code].ACvalue;
#else
				if (IsInHuffmanCodes(code, k,  HTAC->m_numBlocks, HTAC->m_blocks, &decodedValue))

//				if (IsInHuffmanCodes(code, k,  c->m_acTable->m_numBlocks, c->m_acTable->m_blocks, &decodedValue))
				{
#endif
				// Skip over k bits, since we found the huffman value
#if (STREAM_POINTER != 0)
				SkipNBits(&jdata->m_stream, k);
#else
				/********************** FillNBits ***************************/\
					AC_SkipNBits_FillNBits_loop:
					FillNBits(hufUnit,k);
//
				hufUnit->g_nbits_in_reservoir -= (k);
				hufUnit->g_reservoir &= ((1U<<hufUnit->g_nbits_in_reservoir)-1);
//				SkipNBits(jdata->m_stream, k);
#endif


				// Our decoded value is broken down into 2 parts, repeating RLE, and then
				// the number of bits that make up the actual value next
				int valCode = decodedValue;

				unsigned char size_val = valCode&0xF;	// Number of bits for our data
				unsigned char count_0  = valCode>>4;	// Number RunLengthZeros

				if (size_val==0)
				{// RLE
					if (count_0==0)EOB_found=true;	// EOB found, go out
					else if (count_0==0xF) nr+=16;  // skip 16 zeros
				}
				else
				{
					nr+=count_0; //skip count_0 zeroes

					if (nr > 63)
					{
						printf("-|- ##ERROR## Huffman Decoding\n");
//						dprintf("-|- ##ERROR## Huffman Decoding\n");
						///DBG_HALT;
					}

					//if ( jdata->m_stream[0]==0xff && jdata->m_stream[1]!=0x00)
					//{
					//	jdata->m_stream += 2;
					//	g_reservoir = 0;
					//	g_nbits_in_reservoir = 0;
					//}

#if (STREAM_POINTER != 0)
					short data = GetNBits(&jdata->m_stream, size_val );
#else
//					printf("size_val%d\n\n",size_val);
					/********************** FillNBits ***************************/\
						AC_GetNBits_FillNBits_loop:
						FillNBits(hufUnit,size_val);
//
					short data = ((hufUnit->g_reservoir)>>(hufUnit->g_nbits_in_reservoir-(size_val)));

					hufUnit->g_nbits_in_reservoir -= (size_val);
					hufUnit->g_reservoir &= ((1U<<hufUnit->g_nbits_in_reservoir)-1);
//					short data = GetNBits(jdata->m_stream, size_val );
#endif

					data = DetermineSign(data, size_val);
#if HUF_DCT_SPEEDUP == 0
					DCT_tcoeff[nr++]=data;
#else
					DCT[nr++]=data;
#endif
				}
#if HUF_DEBUG == 1
				printf("%d\tbits wide AC Huf Code\n",k);
#endif
//#if HUF_DEL_COUNT == 1
//
//#endif
				break;
			}
		}
#if HUF_DEL_COUNT == 1
//	huffmanACDelay += (k-2) * 162 * ISINHUF_SEARCH_FACTOR + \
//						(huffmanIsInHufCodesdelay * ISINHUF_SEARCH_FACTOR)  + HUF_AC_SEARCH_LATENCY;
	huffmanACDelay += (huffmanIsInHufCodesdelay * ISINHUF_SEARCH_FACTOR + ISINHUF_INIT)  + HUF_AC_SEARCH_LATENCY;
	huffmanIsInHufCodesdelay = 0;
	huffmanAClengthsTotal += k;
#endif
//			if (k>16)
//			{
//				nr++;
//			}
	}
#if HUF_DEBUG == 1
	printf("%d\tnon-zero AC coefficients\n\n",nr);
	huffmanNonZeroCoef += nr;
#if HUF_LUT_SPEEDUP == 0
	printf("%d\t AC m_numBlocks\n",HTAC->m_numBlocks);
#endif
#endif
#if HUF_DEL_COUNT == 1
//	huffmanACDelay = (k-1) * 162 * ISINHUF_SEARCH_DELAY + (nr * ISINHUF_SEARCH_DELAY) + HUF_AC_SEARCH_DELAY;
	huffmanTotalDelay += HUF_INIT_DELAY + huffmanDCDelay + huffmanACDelay;
//	huffmanTotalDelay += 64 + 16*(162*2+9) + (16*(162*2+9)+8)*63;
//	huffmanTotalDelay += 64 + 1*(1*2+9) + (1*(1*2+9)+8)*1;
//	huffmanTotalDelay += 1*1*2 + 1*1*2*63;
	printf("%lu DC Delay + %lu AC Delay = %lu Block Delay\n%lu\tTotal Delay\n\n",\
			huffmanDCDelay, huffmanACDelay, huffmanDCDelay+huffmanACDelay, huffmanTotalDelay);
#endif

//	}


	#ifndef JPG_SPEEDUP
	DumpDCTValues(DCT_tcoeff);
	#endif

#if HUF_DCT_SPEEDUP == 0
	// We've decoded a block of data, so copy it across to our buffer
	Huffman_save_DCT_coef_loop:
	for (int j = 0; j < 64; j++)
	{
//		c->m_DCT[j] = DCT_tcoeff[j];
		m_DCT[j] = DCT_tcoeff[j];
	}
#endif
} // end of ProcessHuffmanBlock()

/***************************************************************************/

inline void ConvertYCrCbtoRGB(int y, int cb, int cr,
							  int* r, int* g, int* b)

{
	float red, green, blue;

	red   = y + 1.402f*(cb-128);
	green = y-0.34414f*(cr-128)-0.71414f*(cb-128);
	blue  = y+1.772f*(cr-128);

	*r = (int) Clamp((int)red);
	*g = (int) Clamp((int)green);
	*b = (int) Clamp((int)blue);
}

/***************************************************************************/

//inline void YCrCB_to_RGB24_Block8x8(stJpegData *jdata, int w, int h, int imgx, int imgy, unsigned int imgw, unsigned int imgh)
inline void YCrCB_to_RGB24_Block8x8(stJpegData *jdata, unsigned char w, unsigned char h, int imgx, int imgy, unsigned int imgw, unsigned int imgh)
{
	const unsigned char *Y, *Cb, *Cr;
	unsigned char *pix;

	int r, g, b;

	Y  = jdata->m_Y;
	Cb = jdata->m_Cb;
	Cr = jdata->m_Cr;

//	int olw = 0; // overlap
////	if ( imgx > abs(imgw-8*w) )
////	{
////		olw = (imgw-imgx)*w + 1;
////	}
//
//	int olh = 0; // overlap
////	if ( imgy > abs(imgh-8*h) )
////	{
////		olh = (imgh-imgy)*2 + 1;
////	}
//
////	dprintf("***pix***\n\n");

//	printf("8*h %d\n", 8*h);
//	printf("8*w %d\n", 8*w);
	YCrCB_to_RGB24_loop_vertical:
//	for (unsigned int y=0; y<(8*h - olh); y++)
	for (unsigned int y=0; y<(8*h); y++)
	{
#pragma HLS LOOP_TRIPCOUNT min=16 max=16 avg=16
//#pragma HLS LOOP_TRIPCOUNT min=1 max=2 avg=1
		YCrCB_to_RGB24_loop_horizontal:
//		for (unsigned int x=0; x<(8*w - olw); x++)
		for (unsigned int x=0; x<(8*w); x++)
		{
#pragma HLS LOOP_TRIPCOUNT min=16 max=16 avg=16
//#pragma HLS LOOP_TRIPCOUNT min=8 max=16 avg=16
#if RGB_OVERF_CHECK == 1
			if (x+imgx >= imgw) continue;
			if (y+imgy >= imgh) continue;
#endif
//			int poff = x*3 + jdata->m_width*3*y;


//			pix = &(jdata->m_colourspace[poff]);

//			printf("poff %d\n", poff);

//			printf("m_colourspace %d\n", pix[13]);

//			printf("m_colourspace %d\n", &jdata->m_colourspace[poff]);
//			printf("m_colourspace SIZE %d\n", sizeof(jdata->m_colourspace));


//			int yoff = x + y*8; //my wrong
//			int yoff = (int)(x*(1.0f/w)) + (int)(y*(1.0f/h));// my wrong
			int yoff = x + y*(w*8);
			int coff = (int)(x*(1.0f/w)) + (int)(y*(1.0f/h))*8;

			int yc =  Y[yoff];
			int cb = Cb[coff];
			int cr = Cr[coff];

			ConvertYCrCbtoRGB(yc,cr,cb,&r,&g,&b);

//			int rgb_index = imgx*3 + (imgy *jdata->m_width*3);
//
//			jdata->m_rgb[rgb_index + x*3 + jdata->m_width*3*y] = Clamp(r);
//			jdata->m_rgb[rgb_index + x*3 + jdata->m_width*3*y + 1] = Clamp(g);
//			jdata->m_rgb[rgb_index + x*3 + jdata->m_width*3*y + 2] = Clamp(b);
			int rgb_index = imgx*3 + (imgy *imgw*3);

			jdata->m_rgb[rgb_index + x*3 + imgw*3*y] = Clamp(r);
			jdata->m_rgb[rgb_index + x*3 + imgw*3*y + 1] = Clamp(g);
			jdata->m_rgb[rgb_index + x*3 + imgw*3*y + 2] = Clamp(b);

//			jdata->m_colourspace[x*3 + jdata->m_width*3*y] = Clamp(r);
//			jdata->m_colourspace[x*3 + jdata->m_width*3*y + 1] = Clamp(g);
//			jdata->m_colourspace[x*3 + jdata->m_width*3*y + 2] = Clamp(b);

			//			jdata->m_colourspace[x*3 + jdata->m_width*3*y] = Clamp(r);
			//			jdata->m_colourspace[x*3 + jdata->m_width*3*y + 1] = Clamp(g);
			//			jdata->m_colourspace[x*3 + jdata->m_width*3*y + 1] = Clamp(b);
//			pix[0] = Clamp(r);
//			pix[1] = Clamp(g);
//			pix[2] = Clamp(b);

//			dprintf("-[%d][%d][%d]-\t", poff, yoff, coff);
		}
//		dprintf("\n");
	}
//	dprintf("\n\n");
}

/***************************************************************************/
//
//  Decoding
//  .-------.
//  | 1 | 2 |
//  |---+---|
//  | 3 | 4 |
//  `-------'
//
/***************************************************************************/
//inline void DecodeMCU(stJpegData *jdata, int w, int h)
inline void DecodeMCU(stJpegData *jdata, unsigned char hFactor, unsigned char vFactor)
{
//	int w = jdata->m_component_info[cY].m_hFactor;
//	int h = jdata->m_component_info[cY].m_vFactor;
	// Y
	DecodeMCU_loop_vertical:
//	for (int y=0; y<h; y++)
//	for (int y=0; y<jdata->m_component_info[cY].m_vFactor; y++)
	for (unsigned char y=0; y<vFactor; y++)
	{
#pragma HLS LOOP_TRIPCOUNT min=2 max=2 avg=2
//#pragma HLS LOOP_TRIPCOUNT min=1 max=2 avg=2
		DecodeMCU_loop_horizontal:
//		for (int x=0; x<w; x++)
//		for (int x=0; x<jdata->m_component_info[cY].m_hFactor; x++)
		for (unsigned char x=0; x<hFactor; x++)
		{
#pragma HLS LOOP_TRIPCOUNT min=2 max=2 avg=2
//#pragma HLS LOOP_TRIPCOUNT min=1 max=2 avg=2
//			int stride = 8;//w*8;
//			int offset = x + y*64;//x*8 + y*64*w;
//			int stride = w*8;
			unsigned char stride = hFactor*8;
//			int offset = x*8 + y*64*w;
			unsigned int offset = x*8 + y*64*hFactor;

//			ProcessHuffmanBlock(jdata, cY);
//			ProcessHuffmanBlock(&jdata->m_component_info[cY], &jdata->m_Huffman, cY);
#if HUF_LUT_SPEEDUP == 1
			ProcessHuffmanBlock(jdata->m_component_info[cY].m_DCT, &jdata->m_component_info[cY].m_previousDC, \
					jdata->m_Huffman.m_hufLUTs[jdata->m_component_info[cY].dcTable_index].m_hufLUT, \
					jdata->m_Huffman.m_hufLUTs[jdata->m_component_info[cY].acTable_index].m_hufLUT, &jdata->m_Huffman);
#else
			ProcessHuffmanBlock(jdata->m_component_info[cY].m_DCT, &jdata->m_component_info[cY].m_previousDC, &jdata->m_Huffman.m_HTDC[jdata->m_component_info[cY].dcTable_index], &jdata->m_Huffman.m_HTAC[jdata->m_component_info[cY].acTable_index], &jdata->m_Huffman);
#endif
//			ProcessHuffmanBlock(jdata->m_component_info[cY], jdata->m_Huffman[cY], m_HTDC[cY], m_HTAC[cY] cY);

			///DBG_ASSERT(cY>0 && cY<COMPONENTS);
			///DBG_ASSERT(offset>=0 && offset<64*4);
//			DecodeSingleBlock(&jdata->m_component_info[cY], &jdata->m_Y[offset], stride);


			DecodeSingleBlock(jdata->m_component_info[cY].m_DCT, jdata->m_component_info[cY].m_qTable, &jdata->m_Y[offset], stride);
		}
	}

	// Cb
//	ProcessHuffmanBlock(jdata, cCb);
//	ProcessHuffmanBlock(&jdata->m_component_info[cCb], &jdata->m_Huffman, cCb);

//	DecodeSingleBlock(&jdata->m_component_info[cCb], jdata->m_Cb, 8);
#if HUF_LUT_SPEEDUP == 1
			ProcessHuffmanBlock(jdata->m_component_info[cCb].m_DCT, &jdata->m_component_info[cCb].m_previousDC, \
					jdata->m_Huffman.m_hufLUTs[jdata->m_component_info[cCb].dcTable_index].m_hufLUT, \
					jdata->m_Huffman.m_hufLUTs[jdata->m_component_info[cCb].acTable_index].m_hufLUT, &jdata->m_Huffman);
#else
			ProcessHuffmanBlock(jdata->m_component_info[cCb].m_DCT, &jdata->m_component_info[cCb].m_previousDC, &jdata->m_Huffman.m_HTDC[jdata->m_component_info[cCb].dcTable_index], &jdata->m_Huffman.m_HTAC[jdata->m_component_info[cCb].acTable_index], &jdata->m_Huffman);
#endif
	DecodeSingleBlock(jdata->m_component_info[cCb].m_DCT, jdata->m_component_info[cCr].m_qTable, jdata->m_Cb, 8);	// , stride = 8
	// Cr

//	ProcessHuffmanBlock(jdata, cCr);
//	ProcessHuffmanBlock(&jdata->m_component_info[cCr], &jdata->m_Huffman, cCr);

//	DecodeSingleBlock(&jdata->m_component_info[cCr], jdata->m_Cr, 8);
#if HUF_LUT_SPEEDUP == 1
			ProcessHuffmanBlock(jdata->m_component_info[cCr].m_DCT, &jdata->m_component_info[cCr].m_previousDC, \
					jdata->m_Huffman.m_hufLUTs[jdata->m_component_info[cCr].dcTable_index].m_hufLUT, \
					jdata->m_Huffman.m_hufLUTs[jdata->m_component_info[cCr].acTable_index].m_hufLUT, &jdata->m_Huffman);
#else
			ProcessHuffmanBlock(jdata->m_component_info[cCr].m_DCT, &jdata->m_component_info[cCr].m_previousDC, &jdata->m_Huffman.m_HTDC[jdata->m_component_info[cCr].dcTable_index], &jdata->m_Huffman.m_HTAC[jdata->m_component_info[cCr].acTable_index], &jdata->m_Huffman);
#endif
	DecodeSingleBlock(jdata->m_component_info[cCr].m_DCT, jdata->m_component_info[cCr].m_qTable, jdata->m_Cr, 8);  // ,stride = 8
}

#if HUF_LUT_SPEEDUP == 1

#if HUF_LUT_NEW == 0
void BuildHuffmanLUT(stJpegData* jdata)
//void BuildHuffmanLUT(stJpegData* jdata, int w, int h)
//void BuildHuffmanLUT(stHuffmanTable* HufTable, int w, int h)
{
		BuildHuffmanLUT_components_loop:
		for (uint8_t comp=0; comp<COMPONENTS; comp++)
		{
			BuildHuffmanLUT_hufVal_key_loop:
			for (uint32_t key_indx=0; key_indx<65536; key_indx++)
			{
#if LUT_SPEEDUP == 0
				jdata->m_HTDC[jdata->m_component_info[comp].dcTable_index].m_hufVal_key[key_indx] = -1;
				jdata->m_HTAC[jdata->m_component_info[comp].acTable_index].m_hufVal_key[key_indx] = -1;
#elif LUT_SPEEDUP == 1
//				jdata->m_Huffman.m_HTDC[jdata->m_component_info[comp].dcTable_index].m_hufVal_key[key_indx] = -1;
//				jdata->m_Huffman.m_HTAC[jdata->m_component_info[comp].acTable_index].m_hufVal_key[key_indx] = -1;
				jdata->m_Huffman.m_hufLUTs[comp].m_hufLUT[key_indx].DClength = 20;
				jdata->m_Huffman.m_hufLUTs[comp].m_hufLUT[key_indx].AClength = 20;
#else
				jdata->m_HTDC[jdata->m_component_info[comp].dcTable_index].m_hufVal_key[key_indx].length = 20;
				jdata->m_HTAC[jdata->m_component_info[comp].acTable_index].m_hufVal_key[key_indx].length = 20;
#endif
			}
//			// save the first huffman code
//			jdata->m_HTDC[jdata->m_component_info[comp].dcTable_index].m_hufVal_key[  jdata->m_HTDC[ jdata->m_component_info[comp].dcTable_index ].m_blocks[0].code  ] = 0;
//			int key_indx_start = 0;
//			int key_indx_end = 0;

			BuildHuffmanLUT_dcBlock_loop:
			for (int dc_block_indx=0; dc_block_indx<jdata->m_Huffman.m_HTDC[jdata->m_component_info[comp].dcTable_index].m_numBlocks; dc_block_indx++)
			{
#pragma HLS LOOP_TRIPCOUNT min=162 max=162 avg=162
//#pragma HLS LOOP_TRIPCOUNT min=30 max=1024 avg=256
//				// fill the LUT with the same huffman code till next code
//				for (int key_indx=key_indx_start; key_indx<jdata->m_HTDC[ jdata->m_component_info[comp].dcTable_index ].m_blocks[0].code; key_indx++)
#if LUT_SPEEDUP == 0
				jdata->m_HTDC[jdata->m_component_info[comp].dcTable_index].m_hufVal_key[  jdata->m_HTDC[ jdata->m_component_info[comp].dcTable_index ].m_blocks[dc_block_indx].code  ] = dc_block_indx;
#elif LUT_SPEEDUP == 1
				jdata->m_Huffman.m_hufLUTs[comp].m_hufLUT[jdata->m_Huffman.m_HTDC[ jdata->m_component_info[comp].dcTable_index ].m_blocks[dc_block_indx].code ].DCvalue
						= jdata->m_Huffman.m_HTDC[ jdata->m_component_info[comp].dcTable_index ].m_blocks[dc_block_indx].value;
				jdata->m_Huffman.m_hufLUTs[comp].m_hufLUT[jdata->m_Huffman.m_HTDC[ jdata->m_component_info[comp].dcTable_index ].m_blocks[dc_block_indx].code ].DClength
						= jdata->m_Huffman.m_HTDC[ jdata->m_component_info[comp].dcTable_index ].m_blocks[dc_block_indx].length;
//				jdata->m_Huffman.m_HTDC[jdata->m_component_info[comp].dcTable_index].m_hufVal_key[  jdata->m_Huffman.m_HTDC[ jdata->m_component_info[comp].dcTable_index ].m_blocks[dc_block_indx].code  ]
//				                                                                          = jdata->m_Huffman.m_HTDC[ jdata->m_component_info[comp].dcTable_index ].m_blocks[dc_block_indx].value;
////				jdata->m_HTDC[jdata->m_component_info[comp].dcTable_index].m_hufVal_key[  jdata->m_HTDC[ jdata->m_component_info[comp].dcTable_index ].m_blocks[dc_block_indx].code  ].length
////				                                                                          = jdata->m_HTDC[ jdata->m_component_info[comp].dcTable_index ].m_blocks[dc_block_indx].length;
#else
				jdata->m_HTDC[jdata->m_component_info[comp].dcTable_index].m_hufVal_key[  jdata->m_HTDC[ jdata->m_component_info[comp].dcTable_index ].m_blocks[dc_block_indx].code  ].value
				                                                                          = jdata->m_HTDC[ jdata->m_component_info[comp].dcTable_index ].m_blocks[dc_block_indx].value;
				jdata->m_HTDC[jdata->m_component_info[comp].dcTable_index].m_hufVal_key[  jdata->m_HTDC[ jdata->m_component_info[comp].dcTable_index ].m_blocks[dc_block_indx].code  ].length
				                                                                          = jdata->m_HTDC[ jdata->m_component_info[comp].dcTable_index ].m_blocks[dc_block_indx].length;

#endif
				printf("DC Huffman Table size: %d\n",jdata->m_Huffman.m_HTDC[jdata->m_component_info[comp].dcTable_index].m_numBlocks);
				printf("DC huffman code: 0x%x\n", jdata->m_Huffman.m_HTDC[ jdata->m_component_info[comp].dcTable_index ].m_blocks[dc_block_indx].code);
				//				printf(" m_hufVal_key [ %d ] = %d\nlength: %d\nvalue: %d\n", jdata->m_HTAC[jdata->m_component_info[comp].acTable_index].m_blocks[block_indx].code , block_indx, jdata->m_HTAC[jdata->m_component_info[comp].acTable_index].m_blocks[block_indx].length, block_indx, jdata->m_HTAC[jdata->m_component_info[comp].acTable_index].m_blocks[block_indx].value);
			}
			printf("\n\n");
			BuildHuffmanLUT_acBlock_loop:
			for (uint8_t ac_block_indx=0; ac_block_indx<jdata->m_Huffman.m_HTAC[jdata->m_component_info[comp].acTable_index].m_numBlocks; ac_block_indx++)
			{
#pragma HLS LOOP_TRIPCOUNT min=162 max=162 avg=162
//#pragma HLS LOOP_TRIPCOUNT min=30 max=1024 avg=256
#if LUT_SPEEDUP == 0
				jdata->m_HTAC[jdata->m_component_info[comp].acTable_index].m_hufVal_key[  jdata->m_HTAC[ jdata->m_component_info[comp].acTable_index ].m_blocks[ac_block_indx].code  ] = ac_block_indx;
#elif LUT_SPEEDUP == 1
				jdata->m_Huffman.m_hufLUTs[comp].m_hufLUT[jdata->m_Huffman.m_HTAC[ jdata->m_component_info[comp].acTable_index ].m_blocks[ac_block_indx].code ].ACvalue
										= jdata->m_Huffman.m_HTAC[ jdata->m_component_info[comp].acTable_index ].m_blocks[ac_block_indx].value;
				jdata->m_Huffman.m_hufLUTs[comp].m_hufLUT[jdata->m_Huffman.m_HTAC[ jdata->m_component_info[comp].acTable_index ].m_blocks[ac_block_indx].code ].AClength
										= jdata->m_Huffman.m_HTAC[ jdata->m_component_info[comp].acTable_index ].m_blocks[ac_block_indx].length;
//				jdata->m_Huffman.m_HTAC[jdata->m_component_info[comp].acTable_index].m_hufVal_key[  jdata->m_HTAC[ jdata->m_component_info[comp].acTable_index ].m_blocks[ac_block_indx].code  ]
//				                                                                          = jdata->m_Huffman.m_HTAC[ jdata->m_component_info[comp].acTable_index ].m_blocks[ac_block_indx].value;
//
////				jdata->m_HTAC[jdata->m_component_info[comp].acTable_index].m_hufVal_key[  jdata->m_HTAC[ jdata->m_component_info[comp].acTable_index ].m_blocks[ac_block_indx].code  ].length
////				                                                                          = jdata->m_HTAC[ jdata->m_component_info[comp].acTable_index ].m_blocks[ac_block_indx].length;
#else
				jdata->m_HTAC[jdata->m_component_info[comp].acTable_index].m_hufVal_key[  jdata->m_HTAC[ jdata->m_component_info[comp].acTable_index ].m_blocks[ac_block_indx].code  ].value
				                                                                          = jdata->m_HTAC[ jdata->m_component_info[comp].acTable_index ].m_blocks[ac_block_indx].value;

				jdata->m_HTAC[jdata->m_component_info[comp].acTable_index].m_hufVal_key[  jdata->m_HTAC[ jdata->m_component_info[comp].acTable_index ].m_blocks[ac_block_indx].code  ].length
				                                                                          = jdata->m_HTAC[ jdata->m_component_info[comp].acTable_index ].m_blocks[ac_block_indx].length;

#endif
				printf("AC Huffman Table size: %d\n",jdata->m_Huffman.m_HTAC[jdata->m_component_info[comp].acTable_index].m_numBlocks);
				printf("AC huffman code: 0x%x\n", jdata->m_Huffman.m_HTAC[ jdata->m_component_info[comp].acTable_index ].m_blocks[ac_block_indx].code);

//				printf(" m_hufVal_key [ %d ] = %d\nlength: %d\nvalue: %d\n", jdata->m_HTAC[jdata->m_component_info[comp].acTable_index].m_blocks[block_indx].code , block_indx, jdata->m_HTAC[jdata->m_component_info[comp].acTable_index].m_blocks[block_indx].length, block_indx, jdata->m_HTAC[jdata->m_component_info[comp].acTable_index].m_blocks[block_indx].value);
			}
		}
//jdata->m_Huffman.m_hufLUTs[0].m_hufLUT[0].DCvalue = 2;
//jdata->m_Huffman.m_hufLUTs[0].m_hufLUT[0].ACvalue = 2;
//jdata->m_Huffman.m_hufLUTs[1].m_hufLUT[0].DCvalue = 2;
//jdata->m_Huffman.m_hufLUTs[1].m_hufLUT[0].ACvalue = 2;
//jdata->m_Huffman.m_hufLUTs[2].m_hufLUT[0].DCvalue = 2;
//jdata->m_Huffman.m_hufLUTs[3].m_hufLUT[0].ACvalue = 2;
}
#else
/***************************************************************************/
//void BuildHuffmanLUT(stJpegData* jdata, int w, int h)
void BuildHuffmanLUT(stJpegData* jdata)
//void BuildHuffmanLUT(stHufLUT *hufLUT, stHuffmanTable *HTDC, stHuffmanTable *HTAC)
//void BuildHuffmanLUT(stHuffmanTable* HufTable, int w, int h)
{
		BuildHuffmanLUT_components_loop:
		for (uint8_t comp=0; comp<COMPONENTS; comp++)
		{
			BuildHuffmanLUT_hufVal_key_loop:
			for (uint32_t key_indx=0; key_indx<65536; key_indx++)
			{
#if LUT_SPEEDUP == 0
//				jdata->m_Huffman.m_HTDC[jdata->m_component_info[comp].dcTable_index].m_hufVal_key[key_indx] = -1;
//				jdata->m_Huffman.m_HTAC[jdata->m_component_info[comp].acTable_index].m_hufVal_key[key_indx] = -1;
#elif LUT_SPEEDUP == 1
//				jdata->m_Huffman.m_HTDC[jdata->m_component_info[comp].dcTable_index].m_hufVal_key[key_indx] = -1;
//				jdata->m_Huffman.m_HTAC[jdata->m_component_info[comp].acTable_index].m_hufVal_key[key_indx] = -1;
				jdata->m_Huffman.m_hufLUTs[comp].m_hufLUT[key_indx].DClength = 20;
				jdata->m_Huffman.m_hufLUTs[comp].m_hufLUT[key_indx].AClength = 20;
#else
//				jdata->m_Huffman.m_HTDC[jdata->m_component_info[comp].dcTable_index].m_hufVal_key[key_indx].length = 20;
//				jdata->m_Huffman.m_HTAC[jdata->m_component_info[comp].acTable_index].m_hufVal_key[key_indx].length = 20;
#endif
			}
//			// save the first huffman code
//			jdata->m_HTDC[jdata->m_component_info[comp].dcTable_index].m_hufVal_key[  jdata->m_HTDC[ jdata->m_component_info[comp].dcTable_index ].m_blocks[0].code  ] = 0;
//			int key_indx_start = 0;
//			int key_indx_end = 0;

			BuildHuffmanLUT_loop:
			for ( uint8_t block_cnt=0; block_cnt<162; block_cnt++)
			{
				if (block_cnt < jdata->m_Huffman.m_HTDC[jdata->m_component_info[comp].dcTable_index].m_numBlocks)
//				if (block_cnt < dc_nof_blocks)
				{
					jdata->m_Huffman.m_hufLUTs[comp].m_hufLUT[ jdata->m_Huffman.m_HTDC[jdata->m_component_info[comp].dcTable_index].m_blocks[block_cnt].code ].DCvalue = jdata->m_Huffman.m_HTDC[jdata->m_component_info[comp].dcTable_index].m_blocks[block_cnt].value;
					jdata->m_Huffman.m_hufLUTs[comp].m_hufLUT[ jdata->m_Huffman.m_HTDC[jdata->m_component_info[comp].dcTable_index].m_blocks[block_cnt].code ].DClength = jdata->m_Huffman.m_HTDC[jdata->m_component_info[comp].dcTable_index].m_blocks[block_cnt].length;
				}
				if (block_cnt < jdata->m_Huffman.m_HTAC[jdata->m_component_info[comp].acTable_index].m_numBlocks)
//				if (block_cnt < ac_nof_blocks)
				{
					jdata->m_Huffman.m_hufLUTs[comp].m_hufLUT[ jdata->m_Huffman.m_HTAC[jdata->m_component_info[comp].acTable_index].m_blocks[block_cnt].code ].ACvalue = jdata->m_Huffman.m_HTAC[jdata->m_component_info[comp].acTable_index].m_blocks[block_cnt].value;
					jdata->m_Huffman.m_hufLUTs[comp].m_hufLUT[ jdata->m_Huffman.m_HTAC[jdata->m_component_info[comp].acTable_index].m_blocks[block_cnt].code ].AClength = jdata->m_Huffman.m_HTAC[jdata->m_component_info[comp].acTable_index].m_blocks[block_cnt].length;
				}
			}
		}

//			BuildHuffmanLUT_dcBlock_loop:
//			for (int dc_block_indx=0; dc_block_indx<jdata->m_Huffman.m_HTDC[jdata->m_component_info[comp].dcTable_index].m_numBlocks; dc_block_indx++)
//			{
//#pragma HLS LOOP_TRIPCOUNT min=162 max=162 avg=162
////#pragma HLS LOOP_TRIPCOUNT min=30 max=1024 avg=256
////				// fill the LUT with the same huffman code till next code
////				for (int key_indx=key_indx_start; key_indx<jdata->m_HTDC[ jdata->m_component_info[comp].dcTable_index ].m_blocks[0].code; key_indx++)
//#if LUT_SPEEDUP == 0
////				jdata->m_Huffman.m_HTDC[jdata->m_component_info[comp].dcTable_index].m_hufVal_key[  jdata->m_Huffman.m_HTDC[ jdata->m_component_info[comp].dcTable_index ].m_blocks[dc_block_indx].code  ] = dc_block_indx;
//#elif LUT_SPEEDUP == 1
////				jdata->m_Huffman.m_HTDC[jdata->m_component_info[comp].dcTable_index].m_hufVal_key[  jdata->m_Huffman.m_HTDC[ jdata->m_component_info[comp].dcTable_index ].m_blocks[dc_block_indx].code  ]
////				                                                                          = jdata->m_Huffman.m_HTDC[ jdata->m_component_info[comp].dcTable_index ].m_blocks[dc_block_indx].value;
//////				jdata->m_HTDC[jdata->m_component_info[comp].dcTable_index].m_hufVal_key[  jdata->m_Huffman.m_HTDC[ jdata->m_component_info[comp].dcTable_index ].m_blocks[dc_block_indx].code  ].length
////				                                                                          = jdata->m_Huffman.m_HTDC[ jdata->m_component_info[comp].dcTable_index ].m_blocks[dc_block_indx].length;
//#else
//				jdata->m_Huffman.m_Huffman.m_HTDC[jdata->m_component_info[comp].dcTable_index].m_hufVal_key[  jdata->m_HTDC[ jdata->m_component_info[comp].dcTable_index ].m_blocks[dc_block_indx].code  ].value
//				                                                                          = jdata->m_HTDC[ jdata->m_component_info[comp].dcTable_index ].m_blocks[dc_block_indx].value;
//				jdata->m_Huffman.m_HTDC[jdata->m_component_info[comp].dcTable_index].m_hufVal_key[  jdata->m_HTDC[ jdata->m_component_info[comp].dcTable_index ].m_blocks[dc_block_indx].code  ].length
//				                                                                          = jdata->m_HTDC[ jdata->m_component_info[comp].dcTable_index ].m_blocks[dc_block_indx].length;
//
//#endif
//				printf("DC Huffman Table size: %d\n",jdata->m_Huffman.m_HTDC[jdata->m_component_info[comp].dcTable_index].m_numBlocks);
//				printf("DC huffman code: 0x%x\n", jdata->m_Huffman.m_HTDC[ jdata->m_component_info[comp].dcTable_index ].m_blocks[dc_block_indx].code);
//				//				printf(" m_hufVal_key [ %d ] = %d\nlength: %d\nvalue: %d\n", jdata->m_HTAC[jdata->m_component_info[comp].acTable_index].m_blocks[block_indx].code , block_indx, jdata->m_HTAC[jdata->m_component_info[comp].acTable_index].m_blocks[block_indx].length, block_indx, jdata->m_HTAC[jdata->m_component_info[comp].acTable_index].m_blocks[block_indx].value);
//			}
//			printf("\n\n");
//			BuildHuffmanLUT_acBlock_loop:
//			for (uint8_t ac_block_indx=0; ac_block_indx<jdata->m_Huffman.m_HTAC[jdata->m_component_info[comp].acTable_index].m_numBlocks; ac_block_indx++)
//			{
//#pragma HLS LOOP_TRIPCOUNT min=162 max=162 avg=162
////#pragma HLS LOOP_TRIPCOUNT min=30 max=1024 avg=256
//#if LUT_SPEEDUP == 0
//				jdata->m_Huffman.m_HTAC[jdata->m_component_info[comp].acTable_index].m_hufVal_key[  jdata->m_Huffman.m_HTAC[ jdata->m_component_info[comp].acTable_index ].m_blocks[ac_block_indx].code  ] = ac_block_indx;
//#elif LUT_SPEEDUP == 1
//				jdata->m_Huffman.m_HTAC[jdata->m_component_info[comp].acTable_index].m_hufVal_key[  jdata->m_Huffman.m_HTAC[ jdata->m_component_info[comp].acTable_index ].m_blocks[ac_block_indx].code  ]
//				                                                                          = jdata->m_Huffman.m_HTAC[ jdata->m_component_info[comp].acTable_index ].m_blocks[ac_block_indx].value;
//
////				jdata->m_HTAC[jdata->m_component_info[comp].acTable_index].m_hufVal_key[  jdata->m_HTAC[ jdata->m_component_info[comp].acTable_index ].m_blocks[ac_block_indx].code  ].length
////				                                                                          = jdata->m_HTAC[ jdata->m_component_info[comp].acTable_index ].m_blocks[ac_block_indx].length;
//#else
//				jdata->m_Huffman.m_HTAC[jdata->m_component_info[comp].acTable_index].m_hufVal_key[  jdata->m_HTAC[ jdata->m_component_info[comp].acTable_index ].m_blocks[ac_block_indx].code  ].value
//				                                                                          = jdata->m_HTAC[ jdata->m_component_info[comp].acTable_index ].m_blocks[ac_block_indx].value;
//
//				jdata->m_Huffman.m_HTAC[jdata->m_component_info[comp].acTable_index].m_hufVal_key[  jdata->m_HTAC[ jdata->m_component_info[comp].acTable_index ].m_blocks[ac_block_indx].code  ].length
//				                                                                          = jdata->m_HTAC[ jdata->m_component_info[comp].acTable_index ].m_blocks[ac_block_indx].length;
//
//#endif
//				printf("AC Huffman Table size: %d\n",jdata->m_Huffman.m_HTAC[jdata->m_component_info[comp].acTable_index].m_numBlocks);
//				printf("AC huffman code: 0x%x\n", jdata->m_Huffman.m_HTAC[ jdata->m_component_info[comp].acTable_index ].m_blocks[ac_block_indx].code);
//
////				printf(" m_hufVal_key [ %d ] = %d\nlength: %d\nvalue: %d\n", jdata->m_HTAC[jdata->m_component_info[comp].acTable_index].m_blocks[block_indx].code , block_indx, jdata->m_HTAC[jdata->m_component_info[comp].acTable_index].m_blocks[block_indx].length, block_indx, jdata->m_HTAC[jdata->m_component_info[comp].acTable_index].m_blocks[block_indx].value);
//			}
//		}
}
#endif
//void BuildHuffmanHash(stJpegData* jdata, int w, int h)
////void BuildHuffmanLUT(stHuffmanTable* HufTable, int w, int h)
//{
//		BuildHuffmanHash_components_loop:
//		for (int comp=0; comp<COMPONENTS; comp++)
//		{
////#pragma HLS LOOP_TRIPCOUNT min=3 max=3 avg=3
//			BuildHuffmanHash_hufVal_key_loop:
//			for (int key_indx=0; key_indx<65536; key_indx++)
//			{
////#pragma HLS LOOP_TRIPCOUNT min=65536 max=65536 avg=65536
//				jdata->m_HTDC[jdata->m_component_info[comp].dcTable_index].m_hufVal_key[key_indx] = -1;
//				jdata->m_HTAC[jdata->m_component_info[comp].acTable_index].m_hufVal_key[key_indx] = -1;
//			}
//			BuildHuffmanHash_acBlock_loop:
//			for (int ac_block_indx=0; ac_block_indx<jdata->m_HTAC[jdata->m_component_info[comp].acTable_index].m_numBlocks; ac_block_indx++)
//			{
//#pragma HLS LOOP_TRIPCOUNT min=1024 max=1024 avg=1024
////#pragma HLS LOOP_TRIPCOUNT min=30 max=1024 avg=256
//				jdata->m_HTAC[jdata->m_component_info[comp].acTable_index].m_hufVal_key[  jdata->m_HTAC[ jdata->m_component_info[comp].acTable_index ].m_blocks[ac_block_indx].code  ] = ac_block_indx;
////				printf(" m_hufVal_key [ %d ] = %d\nlength: %d\nvalue: %d\n", jdata->m_HTAC[jdata->m_component_info[comp].acTable_index].m_blocks[block_indx].code , block_indx, jdata->m_HTAC[jdata->m_component_info[comp].acTable_index].m_blocks[block_indx].length, block_indx, jdata->m_HTAC[jdata->m_component_info[comp].acTable_index].m_blocks[block_indx].value);
//			}
//			BuildHuffmanHash_dcBlock_loop:
//			for (int dc_block_indx=0; dc_block_indx<jdata->m_HTDC[jdata->m_component_info[comp].dcTable_index].m_numBlocks; dc_block_indx++)
//			{
//#pragma HLS LOOP_TRIPCOUNT min=1024 max=1024 avg=1024
////#pragma HLS LOOP_TRIPCOUNT min=30 max=1024 avg=256
//				jdata->m_HTDC[jdata->m_component_info[comp].dcTable_index].m_hufVal_key[  jdata->m_HTDC[ jdata->m_component_info[comp].dcTable_index ].m_blocks[dc_block_indx].code  ] = dc_block_indx;
////				printf(" m_hufVal_key [ %d ] = %d\nlength: %d\nvalue: %d\n", jdata->m_HTAC[jdata->m_component_info[comp].acTable_index].m_blocks[block_indx].code , block_indx, jdata->m_HTAC[jdata->m_component_info[comp].acTable_index].m_blocks[block_indx].length, block_indx, jdata->m_HTAC[jdata->m_component_info[comp].acTable_index].m_blocks[block_indx].value);
//			}
//		}
//}
#endif




/***************************************************************************/

///inline
//int JpegDecode(stJpegData *jdata)
//int JpegDecodeHW(stJpegData *jdata, int jpeg_img_height, int jpeg_img_width, int hFactor, int vFactor)
int JpegDecodeHW(stJpegData *jdata, unsigned int jpeg_img_height, unsigned int jpeg_img_width, unsigned char hFactor, unsigned char vFactor)
{
#if SET_CONST_VALUES != 0
	jpeg_img_height = 512;
	jpeg_img_width = 512;
	hFactor = 2;
	vFactor = 2;
#endif

#if HUF_LUT_SPEEDUP == 1
BuildHuffmanLUT(jdata);
#endif

//	stGnBits gNbits;
	jdata->m_Huffman.g_reservoir = 0;
	jdata->m_Huffman.g_nbits_in_reservoir = 0;
//	g_reservoir = 0;
//	g_nbits_in_reservoir = 0;

////	int hFactor = jdata->m_component_info[cY].m_hFactor;
//	int hFactor = 2; //jdata->m_component_info[cY].m_hFactor;
//	int vFactor = 2; //jdata->m_component_info[cY].m_vFactor;
////	int vFactor = jdata->m_component_info[cY].m_vFactor;

	printf("cY hFactor %u \nvFactor:%u\n",hFactor, vFactor);
//	printf("cCr hFactor %u \nvFactor:%u\n",hFactor, vFactor);
//	printf("cCb hFactor %u \nvFactor:%u\n",hFactor, vFactor);

//	unsigned int hFactor = jdata->m_component_info[cY].m_hFactor;
//	unsigned int vFactor = jdata->m_component_info[cY].m_vFactor;
	// RGB24:
//	if (jdata->m_rgb == NULL)
//	{
//		int h = jdata->m_height*3;
		int h = jpeg_img_height*3;
		int w = jpeg_img_width*3;
//		int w = jdata->m_width*3;
		int height = h + (8*hFactor) - (h%(8*hFactor));
		int width  = w + (8*vFactor) - (w%(8*vFactor));
//		int height = h + (8*hFactor) - (h%(8*hFactor));
//		int width  = w + (8*vFactor) - (w%(8*vFactor));

//		jdata->m_rgb = new unsigned char[width * height + 100]; // 100 is a safetly

//		jdata->m_rgb = new unsigned char[width * height + 100]; // 100 is a safetly

		 printf("width %u \nheight:%u\n",width, height);

//		#ifndef JPG_SPEEDUP
//		memset(jdata->m_rgb, 0, width*height);
//		#endif
//	}
//	if (jdata->m_rgb == NULL)
//	{
//		int h = jdata->m_height*3;
//		int w = jdata->m_width*3;
//		int height = h + (8*hFactor) - (h%(8*hFactor));
//		int width  = w + (8*vFactor) - (w%(8*vFactor));
//
////		jdata->m_rgb = new unsigned char[width * height + 100]; // 100 is a safetly
//
////		jdata->m_rgb = new unsigned char[width * height + 100]; // 100 is a safetly
//
//		 printf("width %d \nheight:%d\n",width, height);
//
////		#ifndef JPG_SPEEDUP
////		memset(jdata->m_rgb, 0, width*height);
////		#endif
//	}

	jdata->m_component_info[0].m_previousDC = 0;
	jdata->m_component_info[1].m_previousDC = 0;
	jdata->m_component_info[2].m_previousDC = 0;
	jdata->m_component_info[3].m_previousDC = 0;

	int xstride_by_mcu = 8*hFactor;
	int ystride_by_mcu = 8*vFactor;

//	// Don't forget to that block can be either 8 or 16 lines
//	unsigned int bytes_per_blocklines = jpeg_img_width*3 * ystride_by_mcu;
//	unsigned int bytes_per_mcu = 3*xstride_by_mcu;

	decode_macroblock_loop_vertical:
	for (unsigned int y=0 ; y<jpeg_img_height; y+=ystride_by_mcu)
	{
#pragma HLS LOOP_TRIPCOUNT min=32 max=32 avg=32		// 512/16=32
		decode_macroblock_loop_horizontal:
		for (unsigned int x=0; x<jpeg_img_width; x+=xstride_by_mcu)
		{
#pragma HLS LOOP_TRIPCOUNT min=32 max=32 avg=32		// 512/16=32
			// Decode MCU Plane
//			DecodeMCU(jdata, hFactor, vFactor);
			DecodeMCU(jdata, hFactor, vFactor);
			YCrCB_to_RGB24_Block8x8(jdata, hFactor, vFactor, x, y, jpeg_img_width, jpeg_img_height);
//			DecodeMCU(jdata);
//			YCrCB_to_RGB24_Block8x8(jdata, x, y, jdata->m_component_info[cY].m_hFactor, jdata->m_component_info[cY].m_vFactor,jpeg_img_width, jpeg_img_height);
		}
	}
#if HUF_DEL_COUNT == 1
	printf("%.2f average k\n",1.0*huffmanAClengthsTotal/huffmanNonZeroCoef);
	printf("%.2f average j\n",1.0*huffmanIsInHufCodesTotal/(huffmanNonZeroCoef));
	printf("%.2f average non-zero AC coef\n",huffmanNonZeroCoef/(32*32*6.0));
	printf("%.2f average latency of ProcessHuffmanBlock\n\n",huffmanTotalDelay/(32*32*6.0));
#endif
	return 0;
}
//
///***************************************************************************/
//int DecodeJpgFileData(const unsigned char* buf, // Jpg file in memory
//					  const int sizeBuf,		// Size jpg in bytes in memory
////					  unsigned char* rgbpix,	// Output rgb pixels
//					  unsigned char** rgbpix,	// Output rgb pixels
//					  unsigned int* width,		// Output image width
//					  unsigned int* height)		// Output image height
//{
//	// Allocate memory for our decoded jpg structure, all our data will be
//	// decompressed and stored in here for the various stages of our jpeg decoding
//
//	//	stJpegData jdec[JDEC_SIZE];					///	size
//	stJpegData* jdec = new stJpegData();
////	memset(jdec, 0, sizeof(stJpegData));
////	for (int i=0; i<JDEC_SIZE; i++)
////	for (int i=0; i<sizeof(stJpegData); i++)
////	{
//////		jdec->m_rgb[i] = 0;
//////		unsigned char*		m_rgb;				// Final Red Green Blue pixel data
//////	//	unsigned char		m_rgb[64*64];		///
//////		unsigned int		m_width;			// Width of image
//////		unsigned int		m_height;			// Height of image
//////
//////		const unsigned char*m_stream;			// Pointer to the current stream
//////		int					m_restart_interval;
//////
//////		stComponent			m_component_info[COMPONENTS];
//////
//////		float				m_Q_tables[COMPONENTS][64];	// quantization tables
//////		stHuffmanTable		m_HTDC[HUFFMAN_TABLES];		// DC huffman tables
//////		stHuffmanTable		m_HTAC[HUFFMAN_TABLES];		// AC huffman tables
//////
//////		// Temp space used after the IDCT to store each components
//////		unsigned char		m_Y[64*4];
//////		unsigned char		m_Cr[64];
//////		unsigned char		m_Cb[64];
//////
//////		// Internal Pointer use for colorspace conversion, do not modify it !!!
//////		unsigned char *		m_colourspace;
////	}
//
//	 printf("jdec size:%d",sizeof(jdec));
//
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
////	for (int i=0; i<RGB_PIX_SIZE; i++)
////	{
//////		rgbpix[i] = jdec->m_rgb[i];
////		rgbpix[i] = jdec->m_rgb[i];
////	}
//	*rgbpix = jdec->m_rgb;
//
//	// Release the memory for our jpeg decoder structure jdec
//	delete jdec;
//
//	return 1;
//}

/***************************************************************************/

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
////	unsigned char *buf;
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
////	dprintf("File size is: %s\n", lengthOfFile+1);
//
//
//	printf("File size is: %d\n", lengthOfFile+1);
//
////	buf = new unsigned char[lengthOfFile + 4];// +4 is safety padding
//
//	unsigned char buf[BUF_SIZE];			///
////	if (buf == NULL)
////	{
////		dprintf("Not enough memory for loading file\n");
////		///DBG_HALT;
////		return 0;
////	}
//	fread(buf, lengthOfFile, 1, fp);
//	fclose(fp);
//
//
//
////	unsigned char rgbpix[RGB_PIX_SIZE];
//	unsigned char* rgbpix = NULL;
//	unsigned int width  = 0;
//	unsigned int height = 0;
//
////	DecodeJpgFileData(buf, lengthOfFile, rgbpix, &width, &height);
//	DecodeJpgFileData(buf, lengthOfFile, &rgbpix, &width, &height);
//
//// printf("rgbpix size:%d\n",sizeof(rgbpix));
//
//	if (rgbpix==NULL)
//	{
//		dprintf("Failed to decode jpg\n");
//		///DBG_HALT;
//		return 0;
//	}
////
////	// Delete our data we read in from the file
////	delete[] buf;
//
//	// Save it
//	WriteBMP24(szBmpFileOutName, width, height, rgbpix);
//
//	 printf("rgbpix size:%d\n",sizeof(rgbpix));
//	// Since we don't need the pixel information anymore, we must
//	// release this as well
////	delete[] rgbpix;
//
//	return 1;
//}
