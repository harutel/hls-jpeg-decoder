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

///***************************************************************************/
////
//// Load one jpeg image, and decompress it, and save the result.
////
///***************************************************************************/

int main()
{
	char jpg_in[] = "../../../../data/Lenna.jpg";
	char bmp_out[] = "../../../../data/out.bmp";

	ConvertJpgFile(jpg_in, bmp_out);

	return (0);
}


