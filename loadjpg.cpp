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


#define JPG_SPEEDUP

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


/***************************************************************************/

float C(int u)
{
    if (u == 0)
         return (1.0f/sqrtf(2));
    else
         return 1.0f;
}


int IDCT_calc(int x, int y, short int block[8][8])
{
//#pragma HLS ARRAY_PARTITION variable=block complete dim=1
	const float PI = 3.14f;
    float sum=0;

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
         }
    }
    return (int) (0.25 * sum);
}

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

void DeZigZag(short int outBlock[64], short int inBlock[64])
{
	DeZigZag_loop:
	for(int i=0; i<64; i++)
	{
		outBlock[ i ] = inBlock[ZigZagArray[i]];
	}
}

/***************************************************************************/

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

inline void DecodeSingleBlock(short int compDCT[64], float compqTable[64], unsigned char outputBuf[64*4], unsigned char stride)
{
	short int coef[64];

	// Copy our data into the temp array
	Decode_SingleBlock_data_cp_loop:
	for (int i=0; i<64; i++)
	{

		coef[i] = compDCT[i];
	}

	// De-Quantize
	short int data[64];
	DequantizeBlock(data, coef, compqTable);

	// De-Zig-Zag
	short int block[64];
	DeZigZag(block, data);

	// Create an 8x8 array
	short int arrayBlock[8][8];
	TransformArray(arrayBlock, block);

	// Inverse DCT
	short int val[8][8];
	PerformIDCT(val, arrayBlock);

	// Level Shift each element (i.e. add 128), and copy to our
	SingleBlock_Save_vertical_loop:
	for (int y = 0; y < 8; y++)
	{
		SingleBlock_Save_horizontal_loop:
		for (int x=0; x<8; x++)
		{
			outputBuf[x] = Clamp(val[x][y]);
		}

		outputBuf += stride;
	}

	#ifndef JPG_SPEEDUP
	DumpDecodedBlock(val);
	#endif
}

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

void  FillNBits(unsigned char **stream, int nbits_wanted)						\
{
	int lll;			\
	while ((int)g_nbits_in_reservoir<nbits_wanted)
	{
		/**printf("*(*stream) is:%d\n",*(*stream));	*/
		unsigned char c = *(*stream);
		(*stream)+=1;						\
		ll++;									\

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

#endif


#else

inline void FillNBits(const unsigned char** stream, int& nbits_wanted)
{
	while ((int)g_nbits_in_reservoir<nbits_wanted)
	{
		/*const */unsigned char c = *(*stream)++;
		g_reservoir <<= 8;

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

#if HUF_DEL_COUNT == 1
	uint8_t j;
#endif
//	for (int j=0; j<162; j++)
	IsinHuffmanCodes_loop:
#if HUF_DEL_COUNT == 1
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


	// Second, the 63 AC coefficient
	int nr=1;
	bool EOB_found=false;


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

#if HUF_LUT_SPEEDUP == 1
			if (AC_hufLUT[code].AClength == k)
			{
				decodedValue = AC_hufLUT[code].ACvalue;
#else
				if (IsInHuffmanCodes(code, k,  HTAC->m_numBlocks, HTAC->m_blocks, &decodedValue))

				{
#endif
				// Skip over k bits, since we found the huffman value
#if (STREAM_POINTER != 0)
				SkipNBits(&jdata->m_stream, k);
#else
				/********************** FillNBits ***************************/\
					AC_SkipNBits_FillNBits_loop:
					FillNBits(hufUnit,k);

				hufUnit->g_nbits_in_reservoir -= (k);
				hufUnit->g_reservoir &= ((1U<<hufUnit->g_nbits_in_reservoir)-1);

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
					}



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

inline void YCrCB_to_RGB24_Block8x8(stJpegData *jdata, unsigned char w, unsigned char h, int imgx, int imgy, unsigned int imgw, unsigned int imgh)
{
	const unsigned char *Y, *Cb, *Cr;
	unsigned char *pix;

	int r, g, b;

	Y  = jdata->m_Y;
	Cb = jdata->m_Cb;
	Cr = jdata->m_Cr;


	YCrCB_to_RGB24_loop_vertical:
	for (unsigned int y=0; y<(8*h); y++)
	{
#pragma HLS LOOP_TRIPCOUNT min=16 max=16 avg=16
//#pragma HLS LOOP_TRIPCOUNT min=1 max=2 avg=1
		YCrCB_to_RGB24_loop_horizontal:
		for (unsigned int x=0; x<(8*w); x++)
		{
#pragma HLS LOOP_TRIPCOUNT min=16 max=16 avg=16
//#pragma HLS LOOP_TRIPCOUNT min=8 max=16 avg=16
#if RGB_OVERF_CHECK == 1
			if (x+imgx >= imgw) continue;
			if (y+imgy >= imgh) continue;
#endif

			int yoff = x + y*(w*8);
			int coff = (int)(x*(1.0f/w)) + (int)(y*(1.0f/h))*8;

			int yc =  Y[yoff];
			int cb = Cb[coff];
			int cr = Cr[coff];

			ConvertYCrCbtoRGB(yc,cr,cb,&r,&g,&b);


			int rgb_index = imgx*3 + (imgy *imgw*3);

			jdata->m_rgb[rgb_index + x*3 + imgw*3*y] = Clamp(r);
			jdata->m_rgb[rgb_index + x*3 + imgw*3*y + 1] = Clamp(g);
			jdata->m_rgb[rgb_index + x*3 + imgw*3*y + 2] = Clamp(b);


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

inline void DecodeMCU(stJpegData *jdata, unsigned char hFactor, unsigned char vFactor)
{
	// Y
	DecodeMCU_loop_vertical:
	for (unsigned char y=0; y<vFactor; y++)
	{
#pragma HLS LOOP_TRIPCOUNT min=2 max=2 avg=2
//#pragma HLS LOOP_TRIPCOUNT min=1 max=2 avg=2
		DecodeMCU_loop_horizontal:
		for (unsigned char x=0; x<hFactor; x++)
		{
#pragma HLS LOOP_TRIPCOUNT min=2 max=2 avg=2
//#pragma HLS LOOP_TRIPCOUNT min=1 max=2 avg=2
			unsigned char stride = hFactor*8;
			unsigned int offset = x*8 + y*64*hFactor;

#if HUF_LUT_SPEEDUP == 1
			ProcessHuffmanBlock(jdata->m_component_info[cY].m_DCT, &jdata->m_component_info[cY].m_previousDC, \
					jdata->m_Huffman.m_hufLUTs[jdata->m_component_info[cY].dcTable_index].m_hufLUT, \
					jdata->m_Huffman.m_hufLUTs[jdata->m_component_info[cY].acTable_index].m_hufLUT, &jdata->m_Huffman);
#else
			ProcessHuffmanBlock(jdata->m_component_info[cY].m_DCT, &jdata->m_component_info[cY].m_previousDC, &jdata->m_Huffman.m_HTDC[jdata->m_component_info[cY].dcTable_index], &jdata->m_Huffman.m_HTAC[jdata->m_component_info[cY].acTable_index], &jdata->m_Huffman);
#endif


			DecodeSingleBlock(jdata->m_component_info[cY].m_DCT, jdata->m_component_info[cY].m_qTable, &jdata->m_Y[offset], stride);
		}
	}

	// Cb

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

				jdata->m_Huffman.m_hufLUTs[comp].m_hufLUT[key_indx].DClength = 20;
				jdata->m_Huffman.m_hufLUTs[comp].m_hufLUT[key_indx].AClength = 20;
#else
				jdata->m_HTDC[jdata->m_component_info[comp].dcTable_index].m_hufVal_key[key_indx].length = 20;
				jdata->m_HTAC[jdata->m_component_info[comp].acTable_index].m_hufVal_key[key_indx].length = 20;
#endif
			}
//			// save the first huffman code


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
}
#else
/***************************************************************************/
void BuildHuffmanLUT(stJpegData* jdata)
{
		BuildHuffmanLUT_components_loop:
		for (uint8_t comp=0; comp<COMPONENTS; comp++)
		{
			BuildHuffmanLUT_hufVal_key_loop:
			for (uint32_t key_indx=0; key_indx<65536; key_indx++)
			{
#if LUT_SPEEDUP == 0
				jdata->m_Huffman.m_HTDC[jdata->m_component_info[comp].dcTable_index].m_hufVal_key[key_indx] = -1;
				jdata->m_Huffman.m_HTAC[jdata->m_component_info[comp].acTable_index].m_hufVal_key[key_indx] = -1;
#elif LUT_SPEEDUP == 1
				jdata->m_Huffman.m_hufLUTs[comp].m_hufLUT[key_indx].DClength = 20;
				jdata->m_Huffman.m_hufLUTs[comp].m_hufLUT[key_indx].AClength = 20;
#else
				jdata->m_Huffman.m_HTDC[jdata->m_component_info[comp].dcTable_index].m_hufVal_key[key_indx].length = 20;
				jdata->m_Huffman.m_HTAC[jdata->m_component_info[comp].acTable_index].m_hufVal_key[key_indx].length = 20;
#endif
			}
//			// save the first huffman code


			BuildHuffmanLUT_loop:
			for ( uint8_t block_cnt=0; block_cnt<162; block_cnt++)
			{
				if (block_cnt < jdata->m_Huffman.m_HTDC[jdata->m_component_info[comp].dcTable_index].m_numBlocks)

				{
					jdata->m_Huffman.m_hufLUTs[comp].m_hufLUT[ jdata->m_Huffman.m_HTDC[jdata->m_component_info[comp].dcTable_index].m_blocks[block_cnt].code ].DCvalue = jdata->m_Huffman.m_HTDC[jdata->m_component_info[comp].dcTable_index].m_blocks[block_cnt].value;
					jdata->m_Huffman.m_hufLUTs[comp].m_hufLUT[ jdata->m_Huffman.m_HTDC[jdata->m_component_info[comp].dcTable_index].m_blocks[block_cnt].code ].DClength = jdata->m_Huffman.m_HTDC[jdata->m_component_info[comp].dcTable_index].m_blocks[block_cnt].length;
				}
				if (block_cnt < jdata->m_Huffman.m_HTAC[jdata->m_component_info[comp].acTable_index].m_numBlocks)

				{
					jdata->m_Huffman.m_hufLUTs[comp].m_hufLUT[ jdata->m_Huffman.m_HTAC[jdata->m_component_info[comp].acTable_index].m_blocks[block_cnt].code ].ACvalue = jdata->m_Huffman.m_HTAC[jdata->m_component_info[comp].acTable_index].m_blocks[block_cnt].value;
					jdata->m_Huffman.m_hufLUTs[comp].m_hufLUT[ jdata->m_Huffman.m_HTAC[jdata->m_component_info[comp].acTable_index].m_blocks[block_cnt].code ].AClength = jdata->m_Huffman.m_HTAC[jdata->m_component_info[comp].acTable_index].m_blocks[block_cnt].length;
				}
			}
		}


}
#endif

#endif




/***************************************************************************/

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


	jdata->m_Huffman.g_reservoir = 0;
	jdata->m_Huffman.g_nbits_in_reservoir = 0;
	printf("cY hFactor %u \nvFactor:%u\n",hFactor, vFactor);
		int h = jpeg_img_height*3;
		int w = jpeg_img_width*3;
		int height = h + (8*hFactor) - (h%(8*hFactor));
		int width  = w + (8*vFactor) - (w%(8*vFactor));

		 printf("width %u \nheight:%u\n",width, height);


	jdata->m_component_info[0].m_previousDC = 0;
	jdata->m_component_info[1].m_previousDC = 0;
	jdata->m_component_info[2].m_previousDC = 0;
	jdata->m_component_info[3].m_previousDC = 0;

	int xstride_by_mcu = 8*hFactor;
	int ystride_by_mcu = 8*vFactor;

//	// Don't forget to that block can be either 8 or 16 lines

	decode_macroblock_loop_vertical:
	for (unsigned int y=0 ; y<jpeg_img_height; y+=ystride_by_mcu)
	{
#pragma HLS LOOP_TRIPCOUNT min=32 max=32 avg=32		// 512/16=32
		decode_macroblock_loop_horizontal:
		for (unsigned int x=0; x<jpeg_img_width; x+=xstride_by_mcu)
		{
#pragma HLS LOOP_TRIPCOUNT min=32 max=32 avg=32		// 512/16=32
			// Decode MCU Plane

			DecodeMCU(jdata, hFactor, vFactor);
			YCrCB_to_RGB24_Block8x8(jdata, hFactor, vFactor, x, y, jpeg_img_width, jpeg_img_height);
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
