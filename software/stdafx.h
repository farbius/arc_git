/*
 * stdafx.h
 *
 *  Created on: Oct 29, 2019
 *      Author: zello
 */

#ifndef SRC_STDAFX_H_
#define SRC_STDAFX_H_

// types:

typedef unsigned int uint;
typedef unsigned char uchar;


#define FFT_WINDOW 	4096


typedef struct {
	uint 	Nwindow;
	uint 	Width;
	float 	delta1;
	float 	delta2;
	uint 	Ncptrd;
	uint 	MaxTracks;
	uint 	MaxTraces;
	float 	clutterTrh;
	uint 	Nlin;
} st_settings;

typedef struct {
	float 	f 			[FFT_WINDOW];
	float 	Sf 			[FFT_WINDOW];
	float 	Sf_fltd 	[FFT_WINDOW];
	float 	des 		[FFT_WINDOW];
	uint 	Begs 		[FFT_WINDOW];
	uint 	Fins 		[FFT_WINDOW];
} st_fpgadata;

#endif /* SRC_STDAFX_H_ */
