/*
 * Soft_Part_II.h
 *
 *  Created on: Oct 30, 2019
 *      Author: zello
 */

#ifndef SRC_SOFT_PART_II_H_
#define SRC_SOFT_PART_II_H_

#include "Soft_Part_I.h"


typedef struct {
	uint 	Res_i;
	uint 	Num_in_Res;
	float 	Center;
	uchar 	Clutter;
} st_trace;


extern st_trace **Soft_Part_II (st_Results *Results, uint Nstep, uint *Ntraces, float deltas1, float deltas2, uint MaxTracks, uint Ncptrd, float clutterThr, uint Nlin);



#endif /* SRC_SOFT_PART_II_H_ */
