#ifndef OPENJPG_H
#define OPENJPG_H

#include "loadjpg.h"


#define HUFF_VALUES_RANGE		257

typedef float quantTable_t;

struct stImageInfo
{
	unsigned int	m_width;			// Width of image
	unsigned int	m_height;			// Height of image
	unsigned char	m_hFactor[COMPONENTS];
	unsigned char	m_vFactor[COMPONENTS];
};
// Takes the rgb pixel values and creates a bitmap file
inline void WriteBMP24(const char* szBmpFileName, int Width, int Height, unsigned char* RGB);


// Pass in the filename and it creates a bitmap
int ConvertJpgFile(char* szJpgFileInName, char * szBmpFileOutName);

// Pass in the whole jpg file from memory, and it decodes it to RGB pixel data
inline int DecodeJpgFileData(const unsigned char* buf, int sizeBuf, unsigned char** rgbpix, unsigned int* width, unsigned int* height);
// Don't forget to delete[] rgbpix if you use DecodeJpgFileData..




#endif //OPENJPG_H
