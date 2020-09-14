#ifndef SRC_SOFT_PART_III_H_
#define SRC_SOFT_PART_III_H_


#include "stdafx.h"



typedef struct {
	float 	Fmean;
	float 	dFmax;
	uint 	Length;
	float 	Curvature;
	float 	MSE;
	uchar	isBreak;
} st_TracksInfo;

extern st_TracksInfo *Soft_Part_III (st_Results *Results, st_trace **traces, uint Nstep, uint Ntraces, uint *Lengths);


#endif /* SRC_SOFT_PART_III_H_ */
