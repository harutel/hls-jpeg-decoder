/*
 *
 * @ Author: Harutyun Elbakyan
 * @ Date:	31.03.2017
 * @ Description: This file opens jpeg file and fills appropriate buffers for jpeg decompression
 *
 * @ Changelog:
 *
 * 	Version: 1.0
 */


#include <stdio.h>		// sprintf(..), fopen(..)
#include <string.h>		// memset(..)

#include "loadjpg.h"
#include "openjpg.h"


/***************************************************************************/

inline void PrintSOF(const unsigned char *stream)
{
	int width;
	int height;
	int nr_components;
	int precision;

	const char *nr_components_to_string[] =	{	"????",
												"Grayscale",
												"????",
												"YCbCr",
												"CYMK" };

	precision = stream[2];
	height = BYTE_TO_WORD(stream+3);
	width  = BYTE_TO_WORD(stream+5);
	nr_components = stream[7];

	printf("> SOF marker\n");
	printf("Size:%dx%d nr_components:%d (%s)  precision:%d\n",
							width, height,
							nr_components,
							nr_components_to_string[nr_components],
							precision);

//	dprintf("> SOF marker\n");
//	dprintf("Size:%dx%d nr_components:%d (%s)  precision:%d\n",
//							width, height,
//							nr_components,
//							nr_components_to_string[nr_components],
//							precision);
}

void GenHuffCodes( int num_codes, stBlock* arr, unsigned char* huffVal )
{
     int hufcounter = 0;
     int codelengthcounter = 1;


     for(int cc=0; cc< num_codes; cc++)
     {
		 while ( arr[cc].length > codelengthcounter )
		 {
			hufcounter = hufcounter << 1;
			codelengthcounter++;
		 }

		 arr[cc].code = hufcounter;
		 arr[cc].value = huffVal[cc];
		 hufcounter = hufcounter + 1;
	}
}

/***************************************************************************/

// Takes two array of bits, and build the huffman table for size, and code

/***************************************************************************/
inline void BuildHuffmanTable(const unsigned char *bits, const unsigned char *stream, stHuffmanTable *HT, unsigned char* hufVal)
{
//	for (int j=1; j<=16; j++)
//	{
//		HT->m_length[j] = bits[j];
//	}

	// Work out the total number of codes
	int numBlocks = 0;
	for (int i=1; i<=16; i++)
	{
//		numBlocks += HT->m_length[i];
		numBlocks += bits[i];
	}
	HT->m_numBlocks = numBlocks;

	// Fill in the data our our blocks, so we know how many bits each
	// one is
	int c=0;
	for (int i=1; i<=16; i++)
	{
//		for (int j=0; j<HT->m_length[i]; j++)
		for (int j=0; j<bits[i]; j++)
		{
			HT->m_blocks[c].length = i;
			c++;
		}

	}

//	GenHuffCodes(HT->m_numBlocks, HT->m_blocks, HT->m_hufVal);
	GenHuffCodes(HT->m_numBlocks, HT->m_blocks, hufVal);
}

/***************************************************************************/

inline void BuildQuantizationTable(quantTable_t *qtable, const unsigned char *ref_table)
{
	int c = 0;

	for (int i=0; i<8; i++)
	{
		for (int j=0; j<8; j++)
		{
			unsigned char val = ref_table[c];

			qtable[c] = val;
			c++;
		}
	}
}

/***************************************************************************/

inline int ParseDQT(stJpegData *jdata, const unsigned char *stream, quantTable_t Q_tables[COMPONENTS][64])
//inline int ParseDQT(stJpegData *jdata, const unsigned char *stream)
{
	int length, qi;
//	float *table;

	printf("> DQT marker\n");
//	dprintf("> DQT marker\n");
	length = BYTE_TO_WORD(stream) - 2;
	stream += 2;	// Skip length

	while (length>0)
	{
		qi = *stream++;

		int qprecision = qi>>4;	 // upper 4 bits specify the precision
		int qindex     = qi&0xf; // index is lower 4 bits

		if (qprecision)
		{
			// precision in this case is either 0 or 1 and indicates the precision
			// of the quantized values;
			// 8-bit (baseline) for 0 and  up to 16-bit for 1

			printf("Error - 16 bits quantization table is not supported\n");
//			dprintf("Error - 16 bits quantization table is not supported\n");
			///DBG_HALT;
		}

		if (qindex>=4)
		{
			printf("Error - No more 4 quantization table is supported (got %d)\n", qi);
//			dprintf("Error - No more 4 quantization table is supported (got %d)\n", qi);
			///DBG_HALT;
		}

		// The quantization table is the next 64 bytes
////		table = jdata->m_Q_tables[qindex];
//		table = (float*)&Q_tables[qindex];

		// the quantization tables are stored in zigzag format, so we
		// use this functino to read them all in and de-zig zag them
//		BuildQuantizationTable(table, stream);
		BuildQuantizationTable((quantTable_t*)&Q_tables[qindex], stream);
		stream += 64;
		length -= 65;
	}
	return 0;
}

/***************************************************************************/

//inline int ParseSOS(stJpegData *jdata, const unsigned char *stream)
inline int 	ParseSOS(stJpegData *jdata, unsigned char *stream)
//		inline int ParseSOS(stJpegData *jdata, unsigned char stream[STREAM_SIZE])
{
	/*
	SOS		16		0xffd8			Start Of Scan
	Ls		16		2Ns + 6			Scan header length
	Ns		8		1-4				Number of image components
	Csj		8		0-255			Scan Component Selector
	Tdj		4		0-1				DC Coding Table Selector
	Taj		4		0-1				AC Coding Table Selector
	Ss		8		0				Start of spectral selection
	Se		8		63				End of spectral selection
	Ah		4		0				Successive Approximation Bit High
	Ai		4		0				Successive Approximation Bit Low
	*/

	unsigned int nr_components = stream[2];

	printf("> SOS marker\n");
//	dprintf("> SOS marker\n");

	if (nr_components != 3)
	{
		printf("Error - We only support YCbCr image\n");
//		dprintf("Error - We only support YCbCr image\n");
		///DBG_HALT;
	}


	stream += 3;
	for (unsigned int i=0;i<nr_components;i++)
	{
		unsigned int cid   = *stream++;
		unsigned int table = *stream++;

		if ((table&0xf)>=4)
		{
			printf("Error - We do not support more than 2 AC Huffman table\n");
//			dprintf("Error - We do not support more than 2 AC Huffman table\n");
			///DBG_HALT;
		}
		if ((table>>4)>=4)
		{
			printf("Error - We do not support more than 2 DC Huffman table\n");
//			dprintf("Error - We do not support more than 2 DC Huffman table\n");
			///DBG_HALT;
		}
		printf("ComponentId:%d  tableAC:%d tableDC:%d\n", cid, table&0xf, table>>4);
//		dprintf("ComponentId:%d  tableAC:%d tableDC:%d\n", cid, table&0xf, table>>4);

//		printf("acTable\t\t%p\n",jdata->m_HTAC[table&0xf]);
//		printf("acTable\taddress\t%p\n",&jdata->m_HTAC[table&0xf]);
		printf("acTable\t\t%p\n",jdata->m_Huffman.m_HTAC[table&0xf]);
		printf("acTable\taddress\t%p\n",&jdata->m_Huffman.m_HTAC[table&0xf]);

		printf("table\t%d\n",table);
		printf("table&0xf\t%d\n",table&0xf);
		printf("table>>4\t%d\n",table>>4);


//		jdata->m_component_info[cid].m_acTable = &jdata->m_HTAC[table&0xf];

		jdata->m_component_info[cid].acTable_index = table&0xf;
		jdata->m_component_info[cid].dcTable_index = table>>4;


//		for(int cnt=0; cnt<HUFFMAN_TABLES; cnt++)
//		{
////			jdata->m_component_info[cid].m_acTable = 1;//jdata->m_HTAC[table&0xf][cnt];
//		unsigned char	m_length[17];		// 17 values from jpg file,
//											// k =1-16 ; L[k] indicates the number of Huffman codes of length k
//		unsigned char	m_hufVal[257];		// 256 codes read in from the jpeg file
//
//		int				m_numBlocks;
//		stBlock			m_blocks[1024];
//		}
//		jdata->m_component_info[cid].m_dcTable = &jdata->m_HTDC[table>>4];
	}

//	jdata->stream_index = 0;	// this will be the index of m_stream
	jdata->m_Huffman.stream_index = 0;	// this will be the index of m_stream

//	printf("stream size:%d\n\n",nr_components);
////	for(unsigned int i=0; i<STREAM_SIZE; i++)
//	for(unsigned int i=0; i<STREAM_SIZE/10; i++)
//	{
//		jdata->m_stream[i] = stream[i+3];
//	}
//	printf("stream size:%d\n\n",nr_components);
	for(unsigned int i=0;i<STREAM_SIZE;i++)
	{
//		jdata->m_stream[i] = stream[i+3];
		jdata->m_Huffman.m_stream[i] = stream[i+3];
	}
//	jdata->m_stream = stream+3;
//	printf("jdata stream pointer addresss %p\n", jdata->m_stream);
//	printf("jdata stream index %d\n", jdata->stream_index);
	printf("jdata stream pointer addresss %p\n", jdata->m_Huffman.m_stream);
	printf("jdata stream index %d\n", jdata->m_Huffman.stream_index);


	return 0;
}

/***************************************************************************/

//inline int ParseDHT(stJpegData *jdata, const unsigned char *stream)
inline int ParseDHT(stJpegData *jdata, unsigned char *stream)
{
	/*
	u8 0xff
	u8 0xc4 (type of segment)
	u16 be length of segment
	4-bits class (0 is DC, 1 is AC, more on this later)
	4-bits table id
	array of 16 u8 number of elements for each of 16 depths
	array of u8 elements, in order of depth
	*/

	unsigned int count, i;
	unsigned char huff_bits[17];
	int length, index;

	length = BYTE_TO_WORD(stream) - 2;
	stream += 2;	// Skip length

	printf("> DHT marker (length=%d)\n", length);
//	dprintf("> DHT marker (length=%d)\n", length);

	while (length>0)
	{
		index = *stream++;

		// We need to calculate the number of bytes 'vals' will takes
		huff_bits[0] = 0;
		count = 0;
		for (i=1; i<17; i++)
		{
			huff_bits[i] = *stream++;
			count += huff_bits[i];
		}

		if (count > 256)
		{
			printf("Error - No more than 1024 bytes is allowed to describe a huffman table");
//			dprintf("Error - No more than 1024 bytes is allowed to describe a huffman table");
			///DBG_HALT;
		}
		if ( (index &0xf) >= HUFFMAN_TABLES)
		{
			printf("Error - No mode than %d Huffman tables is supported\n", HUFFMAN_TABLES);
//			dprintf("Error - No mode than %d Huffman tables is supported\n", HUFFMAN_TABLES);
			///DBG_HALT;
		}
		printf("Huffman table %s n%d\n", (index&0xf0)?"AC":"DC", index&0xf);
		printf("Length of the table: %d\n", count);
//		dprintf("Huffman table %s n%d\n", (index&0xf0)?"AC":"DC", index&0xf);
//		dprintf("Length of the table: %d\n", count);

		if (index & 0xf0 )
		{
////			unsigned char* huffval = jdata->m_HTAC[index&0xf].m_hufVal;
//			unsigned char* huffval = jdata->m_Huffman.m_HTAC[index&0xf].m_hufVal;
//			unsigned char* huffval ; //= jdata->m_Huffman.m_HTAC[index&0xf].m_hufVal;
			unsigned char huffval[HUFF_VALUES_RANGE] ; //= jdata->m_Huffman.m_HTAC[index&0xf].m_hufVal;
			for (i = 0; i < count; i++)
				huffval[i] = *stream++;

////			BuildHuffmanTable(huff_bits, stream, &jdata->m_HTAC[index&0xf]); // AC
//			BuildHuffmanTable(huff_bits, stream, &jdata->m_Huffman.m_HTAC[index&0xf]); // AC
			BuildHuffmanTable(huff_bits, stream, &jdata->m_Huffman.m_HTAC[index&0xf], huffval); // AC
		}
		else
		{
////			unsigned char* huffval = jdata->m_HTDC[index&0xf].m_hufVal;
//			unsigned char* huffval = jdata->m_Huffman.m_HTDC[index&0xf].m_hufVal;
			unsigned char huffval[HUFF_VALUES_RANGE]; // = jdata->m_Huffman.m_HTDC[index&0xf].m_hufVal;
			for (i = 0; i < count; i++)
				huffval[i] = *stream++;

////			BuildHuffmanTable(huff_bits, stream, &jdata->m_HTDC[index&0xf]); // DC
//			BuildHuffmanTable(huff_bits, stream, &jdata->m_Huffman.m_HTDC[index&0xf]); // DC
			BuildHuffmanTable(huff_bits, stream, &jdata->m_Huffman.m_HTDC[index&0xf], huffval); // DC
		}

		length -= 1;
		length -= 16;
		length -= count;
	}
//	dprintf("< DHT marker\n");
	printf("< DHT marker\n");
	return 0;
}

/***************************************************************************/

//inline
int ParseSOF(stJpegData *jdata, stImageInfo *jinfo, const unsigned char *stream, quantTable_t Q_tables[COMPONENTS][64])
//int ParseSOF(stJpegData *jdata, const unsigned char *stream)
{
	/*
	SOF		16		0xffc0		Start Of Frame
	Lf		16		3Nf+8		Frame header length
	P		8		8			Sample precision
	Y		16		0-65535		Number of lines
	X		16		1-65535		Samples per line
	Nf		8		1-255		Number of image components (e.g. Y, U and V).

	---------Repeats for the number of components (e.g. Nf)-----------------
	Ci		8		0-255		Component identifier
	Hi		4		1-4			Horizontal Sampling Factor
	Vi		4		1-4			Vertical Sampling Factor
	Tqi		8		0-3			Quantization Table Selector.
	*/

	PrintSOF(stream);

	int height = BYTE_TO_WORD(stream+3);
	int width  = BYTE_TO_WORD(stream+5);
	int nr_components = stream[7];



	stream += 8;
	for (int i=0; i<nr_components; i++)
	{
		int cid				= *stream++;
		int sampling_factor = *stream++;
		int Q_table			= *stream++;

		stComponent *c = &jdata->m_component_info[cid];
		jinfo->m_vFactor[cid] = (unsigned char) sampling_factor&0xf;
		jinfo->m_hFactor[cid] = (unsigned char) sampling_factor>>4;
//		c->m_qTable = jdata->m_Q_tables[Q_table];
		for (int cnt=0; cnt<64; cnt++)
		{
//			c->m_qTable[cnt] = jdata->m_Q_tables[Q_table][cnt];
			c->m_qTable[cnt] = Q_tables[Q_table][cnt];
		}

		printf("Component:%d  factor:%ux%u  Quantization table:%d\n",
				cid,
				jinfo->m_vFactor[cid],
				jinfo->m_hFactor[cid],
				Q_table );
//		dprintf("Component:%d  factor:%dx%d  Quantization table:%d\n",
//				cid,
//				c->m_vFactor,
//				c->m_hFactor,
//				Q_table );
	}
//	jdata->m_width = width;
//	jdata->m_height = height;
	jinfo->m_width = width;
	jinfo->m_height = height;

	return 0;
}

/***************************************************************************/

//inline int ParseJFIF(stJpegData *jdata, const unsigned char *stream)
inline int ParseJFIF(stJpegData *jdata, stImageInfo *jinfo, unsigned char *stream)
{
	int chuck_len;
	int marker;
	int sos_marker_found = 0;
	int dht_marker_found = 0;

	quantTable_t quant_tables[COMPONENTS][64];

	// Parse marker
	while (!sos_marker_found)
	{
		if (*stream++ != 0xff)
		{
			goto bogus_jpeg_format;
		}

		// Skip any padding ff byte (this is normal)
		while (*stream == 0xff)
		{
			stream++;
		}

		marker = *stream++;
		chuck_len = BYTE_TO_WORD(stream);

		switch (marker)
		{
			case SOF:
			{
				if (ParseSOF(jdata, jinfo, stream, quant_tables) < 0)
//				if (ParseSOF(jdata, stream) < 0)
					return -1;
			}
			break;

			case DQT:
			{
				if (ParseDQT(jdata, stream, quant_tables) < 0)
//				if (ParseDQT(jdata, stream) < 0)
					return -1;
			}
			break;

			case SOS:
			{
				if (ParseSOS(jdata, stream) < 0)
					return -1;
				sos_marker_found = 1;
			}
			break;

			case DHT:
			{
				if (ParseDHT(jdata, stream) < 0)
					return -1;
				dht_marker_found = 1;
			}
			break;

			// The reason I added these additional skips here, is because for
			// certain jpg compressions, like swf, it splits the encoding
			// and image data with SOI & EOI extra tags, so we need to skip
			// over them here and decode the whole image
			case SOI:
			case EOI:
			{
				chuck_len = 0;
				break;
			}
			break;

			case 0xDD: //DRI: Restart_markers=1;
			{
				jdata->m_Huffman.m_restart_interval = BYTE_TO_WORD(stream);
				printf("DRI - Restart_marker\n");
//				dprintf("DRI - Restart_marker\n");
			}
			break;

			case APP0:
			{
				printf("APP0 Chunk ('txt' information) skipping\n");
//				dprintf("APP0 Chunk ('txt' information) skipping\n");
			}
			break;

			default:
			{
				printf("ERROR> Unknown marker %2.2x\n", marker);
//				dprintf("ERROR> Unknown marker %2.2x\n", marker);
			}
			break;
		}

		stream += chuck_len;
	}

	if (!dht_marker_found)
	{
		printf("ERROR> No Huffman table loaded\n");
//		dprintf("ERROR> No Huffman table loaded\n");
		///DBG_HALT;
	}

	return 0;

	bogus_jpeg_format:
	printf("ERROR> Bogus jpeg format\n");
//	dprintf("ERROR> Bogus jpeg format\n");
	///DBG_HALT;
	return -1;
}

/***************************************************************************/

//inline int JpegParseHeader(stJpegData *jdata, const unsigned char *buf, unsigned int size)
inline int JpegParseHeader(stJpegData *jdata, stImageInfo *jinfo, unsigned char *buf, unsigned int size)
{
	// Identify the file
	if ((buf[0] != 0xFF) || (buf[1] != SOI))
	{
//		dprintf("Not a JPG file ?\n");
		printf("Not a JPG file ?\n");
		///DBG_HALT;
		return -1;
	}

//	const unsigned char* startStream = buf+2;
	unsigned char* startStream = buf+2;
	const int fileSize = size-2;

	printf("-|- File thinks its size is: %d bytes\n", fileSize);
//	dprintf("-|- File thinks its size is: %d bytes\n", fileSize);

	int ret = ParseJFIF(jdata, jinfo, startStream);


//	printf("\n ret: %d\n", ret);

	return ret;
}

/***************************************************************************/
//
// Save a buffer in 24bits Bitmap (.bmp) format
//
/***************************************************************************/
//inline void WriteBMP24(const char* szBmpFileName, int Width, int Height, unsigned char* RGB)
inline void WriteBMP24(const char* szBmpFileName, unsigned int Width, unsigned int Height, unsigned char* RGB)
{
	#pragma pack(1)
	struct stBMFH // BitmapFileHeader & BitmapInfoHeader
	{
		// BitmapFileHeader
		char         bmtype[2];     // 2 bytes - 'B' 'M'
		unsigned int iFileSize;     // 4 bytes
		short int    reserved1;     // 2 bytes
		short int    reserved2;     // 2 bytes
		unsigned int iOffsetBits;   // 4 bytes
		// End of stBMFH structure - size of 14 bytes
		// BitmapInfoHeader
		unsigned int iSizeHeader;    // 4 bytes - 40
		unsigned int iWidth;         // 4 bytes
		unsigned int iHeight;        // 4 bytes
		short int    iPlanes;        // 2 bytes
		short int    iBitCount;      // 2 bytes
		unsigned int Compression;    // 4 bytes
		unsigned int iSizeImage;     // 4 bytes
		unsigned int iXPelsPerMeter; // 4 bytes
		unsigned int iYPelsPerMeter; // 4 bytes
		unsigned int iClrUsed;       // 4 bytes
		unsigned int iClrImportant;  // 4 bytes
		// End of stBMIF structure - size 40 bytes
		// Total size - 54 bytes
	};
	#pragma pack()

	// Round up the width to the nearest DWORD boundary
	int iNumPaddedBytes = 4 - (Width * 3) % 4;
	iNumPaddedBytes = iNumPaddedBytes % 4;

	stBMFH bh;
	memset(&bh, 0, sizeof(bh));
	bh.bmtype[0]='B';
	bh.bmtype[1]='M';
	bh.iFileSize = (Width*Height*3) + (Height*iNumPaddedBytes) + sizeof(bh);
	bh.iOffsetBits = sizeof(stBMFH);
	bh.iSizeHeader = 40;
	bh.iPlanes = 1;
	bh.iWidth = Width;
	bh.iHeight = Height;
	bh.iBitCount = 24;


	char temp[2048]={0};
	sprintf(temp, "%s", szBmpFileName);
	FILE* fp = fopen(temp, "wb");
	///DBG_ASSERT( fp ); // Error creating file - valid file name?  valid location?
	fwrite(&bh, sizeof(bh), 1, fp);
	for (int y=Height-1; y>=0; y--)
	{
		for (int x=0; x<Width; x++)
		{
			int i = (x + (Width)*y) * 3;
			unsigned int rgbpix = (RGB[i]<<16)|(RGB[i+1]<<8)|(RGB[i+2]<<0);
			fwrite(&rgbpix, 3, 1, fp);
		}
		if (iNumPaddedBytes>0)
		{
			unsigned char pad = 0;
			fwrite(&pad, iNumPaddedBytes, 1, fp);
		}
	}
	fclose(fp);
}

/***************************************************************************/
//int DecodeJpgFileData(const unsigned char* buf, // Jpg file in memory
int DecodeJpgFileData(unsigned char* buf, // Jpg file in memory
					  const int sizeBuf,		// Size jpg in bytes in memory
//					  unsigned char* rgbpix,	// Output rgb pixels
					  unsigned char** rgbpix,	// Output rgb pixels
					  unsigned int* width,		// Output image width
					  unsigned int* height)		// Output image height
{
//	// Allocate memory for our decoded jpg structure, all our data will be
//	// decompressed and stored in here for the various stages of our jpeg decoding
//
//	//	stJpegData jdec[JDEC_SIZE];					///	size
//	stJpegData* jdec = new stJpegData();
////	memset(jdec, 0, sizeof(stJpegData));
////	for (int i=0; i<JDEC_SIZE; i++)
////	for (int i=0; i<sizeof(stJpegData); i++)
////	{
////	}
//
//	 printf("jdec size:%d",sizeof(jdec));
//
//	if (jdec == NULL)
//	{
////		dprintf("Not enough memory to alloc the structure need for decompressing\n");
//		printf("Not enough memory to alloc the structure need for decompressing\n");
//		///DBG_HALT;
//		return 0;
//	}
//
//	// Start Parsing.....reading & storing data
//	if (JpegParseHeader(jdec, buf, sizeBuf)<0)
//	{
////		dprintf("ERROR > parsing jpg header\n");
//		printf("ERROR > parsing jpg header\n");
//		///DBG_HALT;
//	}
//
//	// We've read it all in, now start using it, to decompress and create rgb values
//	printf("Decoding JPEG image...\n");
////	dprintf("Decoding JPEG image...\n");
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
//	for (int i=0; i<RGB_PIX_SIZE; i++)
//	{
//		printf("jdec->m_rgb\t\t%d\n",jdec->m_rgb[i]);
//		if (jdec->m_rgb[i] == 0)
//			printf("i\t%d\n",i);
//	}
////	*rgbpix = jdec->m_rgb;
//
//	// Release the memory for our jpeg decoder structure jdec
////	delete jdec;
//
//	return 1;
}

/***************************************************************************/
//
//  Returns the size of the file in bytes
//
/***************************************************************************/
inline int FileSize(FILE *fp)
{
	long pos;
	fseek(fp, 0, SEEK_END);
	pos = ftell(fp);
	fseek(fp, 0, SEEK_SET);
	return pos;
}

/***************************************************************************/
//
// Load one jpeg image, and decompress it, and save the result.
//
/***************************************************************************/
int ConvertJpgFile(char* szJpgFileInName, char * szBmpFileOutName)
{
	FILE *fp;
	unsigned int lengthOfFile;
//	unsigned char *buf;


	// Load the Jpeg into memory
	fp = fopen(szJpgFileInName, "rb");
	if (fp == NULL)
	{
//		dprintf("Cannot open jpg file: %s\n", szJpgFileInName);
		printf("Cannot open jpg file: %s\n", szJpgFileInName);
		///DBG_HALT;
		return 0;
	}

	lengthOfFile = FileSize(fp);
//	dprintf("File size is: %s\n", lengthOfFile+1);


	printf("File size is: %d\n", lengthOfFile+1);

//	buf = new unsigned char[lengthOfFile + 4];// +4 is safety padding

	unsigned char buf[BUF_SIZE];			///
//	if (buf == NULL)
//	{
//		dprintf("Not enough memory for loading file\n");
//		///DBG_HALT;
//		return 0;
//	}
	fread(buf, lengthOfFile, 1, fp);
	fclose(fp);



////	unsigned char rgbpix[RGB_PIX_SIZE];
//
////	unsigned char* rgbpix = NULL;
//	unsigned int width  = 0;
//	unsigned int height = 0;

/*************** JpegDecode **********************/
//	unsigned char* buf, // Jpg file in memory
//						  const int sizeBuf,		// Size jpg in bytes in memory
//	//					  unsigned char* rgbpix,	// Output rgb pixels
//						  unsigned char** rgbpix,	// Output rgb pixels
//						  unsigned int* width,		// Output image width
//						  unsigned int* height)		// Output image height
	// Allocate memory for our decoded jpg structure, all our data will be
	// decompressed and stored in here for the various stages of our jpeg decoding

	//	stJpegData jdec[JDEC_SIZE];					///	size
	stJpegData* jdec = new stJpegData();
//	memset(jdec, 0, sizeof(stJpegData));
//	for (int i=0; i<JDEC_SIZE; i++)
//	for (int i=0; i<sizeof(stJpegData); i++)
//	{
//	}

	 printf("jdec size\t\t\t:%d",sizeof(jdec));

	stImageInfo* jinfo = new stImageInfo();

	if (jdec == NULL)
	{
//		dprintf("Not enough memory to alloc the structure need for decompressing\n");
		printf("Not enough memory to alloc the structure need for decompressing\n");
		///DBG_HALT;
		return 0;
	}

	// Start Parsing.....reading & storing data
	if (JpegParseHeader(jdec, jinfo, buf, lengthOfFile)<0)
	{
//		dprintf("ERROR > parsing jpg header\n");
		printf("ERROR > parsing jpg header\n");
		///DBG_HALT;
	}

	// We've read it all in, now start using it, to decompress and create rgb values
	printf("Decoding JPEG image...\n");
//	dprintf("Decoding JPEG image...\n");

//	JpegDecodeHW(jdec, jinfo->m_height, jinfo->m_width, jdec->m_component_info[cY].m_hFactor, jdec->m_component_info[cY].m_vFactor);
//	JpegDecodeHW(jdec, jinfo->m_height, jinfo->m_width);
	JpegDecodeHW(jdec, jinfo->m_height, jinfo->m_width, jinfo->m_hFactor[cY], jinfo->m_vFactor[cY]);
//	JpegDecode(jdec);
//	// Get the size of the image
//	JpegGetImageSize(jdec, &width, &height);

//	for (int i=0; i<RGB_PIX_SIZE; i++)
//	{
////		rgbpix[i] = jdec->m_rgb[i];
//		rgbpix[i] = jdec->m_rgb[i];
//	}
//	for (int i=0; i<RGB_PIX_SIZE; i++)
//	{
//		printf("jdec->m_rgb\t\t%d\n",jdec->m_rgb[i]);
//		if (jdec->m_rgb[i] == 0)
//			printf("i\t%d\n",i);
//	}


	/*************** JpegDecode **********************/

//
//	DecodeJpgFileData(buf, lengthOfFile, rgbpix, &width, &height);
////	DecodeJpgFileData(buf, lengthOfFile, &rgbpix, &width, &height);

// printf("rgbpix size:%d\n",sizeof(rgbpix));
	if (jdec->m_rgb==NULL)
	{
//		dprintf("Failed to decode jpg\n");
		printf("Failed to decode jpg\n");
		///DBG_HALT;
		return 0;
	}
//
//	// Delete our data we read in from the file
//	delete[] buf;

	// Save it
//	WriteBMP24(szBmpFileOutName, width, height, jdec->m_rgb);
//	WriteBMP24(szBmpFileOutName, jdec->m_width, jdec->m_height, jdec->m_rgb);
	WriteBMP24(szBmpFileOutName, jinfo->m_width, jinfo->m_height, jdec->m_rgb);

	 printf("rgbpix size:%d\n",sizeof(jdec->m_rgb));


		printf("jdata stream pointer addresss %p\n", jdec->m_Huffman.m_stream);

		printf("jdata stream index %d\n", jdec->m_Huffman.stream_index);

		printf("jdata stream pointer addresss %p\n", jdec->m_Huffman.m_stream+jdec->m_Huffman.stream_index);

//	if (rgbpix==NULL)
//	{
////		dprintf("Failed to decode jpg\n");
//		printf("Failed to decode jpg\n");
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
	// Since we don't need the pixel information anymore, we must
	// release this as well
//	delete[] rgbpix;

		delete[] jdec;

	return 1;
}
