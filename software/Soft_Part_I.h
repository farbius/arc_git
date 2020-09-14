#ifndef SRC_SOFT_PART_I_H_
#define SRC_SOFT_PART_I_H_

#include "stdafx.h"

typedef struct {
	uint 	Begs[1024];
	uint 	Fins[1024];
	uchar 	Types[1024];
	float 	Centers[1024];
	uint 	Widths[1024];
	uint 	Num;
	uint 	n_iter;
	char 	isInTrace[1024];
} st_Results;

extern void Soft_Part_I (st_Results *R, float *Sf, uint *Begs, uint *Fins, uint Width);

#endif /* SRC_SOFT_PART_I_H_ */
