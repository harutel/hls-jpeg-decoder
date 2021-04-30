/***************************************************************************/
/*                                                                         */
/*  File: loadjpg.h                                                        */
/*  Author: bkenwright@xbdev.net                                           */
/*  URL: www.xbdev.net                                                     */
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
*/
/***************************************************************************/
#ifndef LOADJPG_H
#define LOADJPG_H

#define SET_CONST_VALUES	(0)

#define RGB_OVERF_CHECK		(1)


#define HUF_DCT_SPEEDUP		(1)



#define PIC					(3)
#define	NEW_HUFFMAN			(0) // it's not working properly


#define FILLNBITS_SPEEDUP	NEW_HUFFMAN
#define HUF_LUT_SPEEDUP		NEW_HUFFMAN
#define HUF_LUT_NEW			(1)
#define LUT_SPEEDUP			(1)

#define	HUF_DEL_COUNT		(1)
#define	HUF_DEBUG			(0 &(HUF_DEL_COUNT))

#define ISINHUF_INIT			(9)
#define	ISINHUF_SEARCH_FACTOR	(2)		// 2 clocks for every read+check
#define	HUF_INIT_DELAY			(64)
#define HUF_DC_SEARCH_LATENCY 	(8)
#define HUF_AC_SEARCH_LATENCY 	(8)

//#define JPG_FILE_SIZE		6//6//1730						// size in KB
//#define IMG_MAX_WIDTH		200 //1280
//#define IMG_MAX_HEIGHT		200//720

#define JPG_FILE_SIZE		105						// size in KB
#define IMG_MAX_WIDTH		512
#define IMG_MAX_HEIGHT		512

//#define JPG_FILE_SIZE		1100//6//1730						// size in KB
//#define IMG_MAX_WIDTH		3264 //1280
//#define IMG_MAX_HEIGHT		2448 //720


#define BUF_SIZE			JPG_FILE_SIZE * 1040 		// something for than actual file size in bytes
#define JDEC_SIZE			110 //105000
//#define RGB_PIX_SIZE		8
#define RGB_PIX_SIZE		3*(IMG_MAX_WIDTH * IMG_MAX_HEIGHT) + 100
#define COLOURSPACE_SIZE	3*RGB_PIX_SIZE
#define STREAM_SIZE			JPG_FILE_SIZE * 1000 		// approx. as actual file size in bytes

#define DQT 	 0xDB	// Define Quantization Table
#define SOF 	 0xC0	// Start of Frame (size information)
#define DHT 	 0xC4	// Huffman Table
#define SOI 	 0xD8	// Start of Image
#define SOS 	 0xDA	// Start of Scan
#define EOI 	 0xD9	// End of Image, or End of File
#define APP0	 0xE0

#define BYTE_TO_WORD(x) (((x)[0]<<8)|(x)[1])


#define COMPONENTS				4
#define HUFFMAN_TABLES			COMPONENTS

#define		cY		1
#define		cCb		2
#define		cCr	 	3
#define		DC		0
#define		AC		1

struct stGnBits
{
	unsigned int g_nbits_in_reservoir;
	unsigned int g_reservoir;
};

struct stBlock
{
     int value;					// Decodes to.
     int length;				// Length in bits.
     unsigned short int code;	// 2 byte code (variable length)
};

/***************************************************************************/


struct stHuffmanTable
{
	int				m_numBlocks;
	stBlock			m_blocks[1024];
};


struct stComponent
{
  float					m_qTable[64];			// quantisation table
////  stHuffmanTable*		m_acTable;
  unsigned int			acTable_index;		// index of HTAC table
  unsigned int 			dcTable_index;		// index of HTDC table

  short int				m_DCT[64];			// DCT coef
  short int				m_previousDC;
};

struct stHufLUT
{
	unsigned char DClength;
	int DCvalue;
	unsigned char AClength;
	int ACvalue;
};

struct stHufLUTs
{
	stHufLUT			m_hufLUT[65536];
};

struct stHuffmanData
{
	stHuffmanTable		m_HTDC[HUFFMAN_TABLES];		// DC huffman tables
	stHuffmanTable		m_HTAC[HUFFMAN_TABLES];		// AC huffman tables
#if HUF_LUT_SPEEDUP == 1
//	unsigned char	m_hufVal_key[65536];
////#elif LUT_SPEEDUP == 1
//	stHufLUT		m_hufVal_key[65536*65536];
//#else
	stHufLUTs			m_hufLUTs[HUFFMAN_TABLES];
#endif
	unsigned char  		m_stream[STREAM_SIZE];			// Pointer to the current stream
	unsigned int		stream_index;
	int 				m_restart_interval;
	unsigned int 		g_nbits_in_reservoir;
	unsigned int 		g_reservoir;

};

struct stJpegData
{
//	unsigned char*		m_rgb;				// Final Red Green Blue pixel data
	unsigned char		m_rgb[RGB_PIX_SIZE];		///

	stComponent			m_component_info[COMPONENTS];

	// Temp space used after the IDCT to store each components
	unsigned char		m_Y[64*4];
	unsigned char		m_Cr[64];
	unsigned char		m_Cb[64];

//	// Internal Pointer use for colorspace conversion, do not modify it !!!
	stHuffmanData		m_Huffman;
};

#ifndef DBG_HALT

	#define DBG_HALT __asm{ int 0x10 }
	#define DBG_ASSERT(exp) {if ( !(exp) ) {DBG_HALT;}}

#endif //DBG_HALT

int JpegDecodeHW(stJpegData *jdata, unsigned int jpeg_img_height, unsigned int jpeg_img_width, unsigned char hFactor, unsigned char vFactor);


void JpegGetImageSize(stJpegData *jdata, unsigned int *width, unsigned int *height);

// Pass in the whole jpg file from memory, and it decodes it to RGB pixel data
inline int DecodeJpgFileData(const unsigned char* buf, int sizeBuf, unsigned char** rgbpix, unsigned int* width, unsigned int* height);
// Don't forget to delete[] rgbpix if you use DecodeJpgFileData..


#endif //LOADJPG_H
