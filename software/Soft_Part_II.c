#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "stdafx.h"
#include "Soft_Part_I.h"
#include "Soft_Part_II.h"

// Local function:
void Soft_Part_II_filtering (st_Results *Results, uint Nstep, float deltas1, float deltas2);
uint GetNextClosest (st_Results *R, uint k, uint *n2, uchar *isClutter, float delta1, float delta2, float ClutterThr);
void SetWholeTraceTo (st_Results *R, st_trace *trace, const uint length, const int num);
float Soft_Part_II_interp1 (st_trace **traces, uint stepnum, uint tracenum, uint Nlin);

// Functions code:

st_trace **Soft_Part_II (st_Results *Results, uint Nstep, uint *Ntraces, float deltas1, float deltas2, uint MaxTracks, uint Ncptrd, float clutterThr, uint Nlin) {

	// 0. variables:
	uint 		i, j, k;
	uint 		i2, k2;
	uint 		i_prev, k_prev;
	st_trace 	**traces;
	st_trace 	*trace;
	//uint 		Ntraces = 0;
	uint 		trace_length;
	uchar 		flag_too_many_traces = 0;
	uchar 		isClutter;
	//uchar 		isClutter_prev;
	uint 		tmp;

	*Ntraces = 0;

	// 1. isInTrace reset:
	for (i = 0; i < Nstep; i++) {
		for (j = 0; j < 1024; j++) Results[i].isInTrace[j] = -1;
	}


	// 2. filtering of a mark:
	Soft_Part_II_filtering (Results, Nstep, deltas1, deltas2);


	// 3. a traces:
	traces = malloc(sizeof(st_trace *) * Nstep);
	for (i = 0; i < Nstep; i++) {
		traces[i] = calloc(sizeof(st_trace), MaxTracks);
	}

	trace = malloc(sizeof(st_trace) * Nstep);

	for (i = 1; i <= Nstep-2; i++) {
		for (k = 1; k <= Results[i-1].Num; k++) {

			if (Results[i-1].isInTrace[k-1] == 0) {

				// cleaning trace:
				for (j = 0; j < Nstep; j++) {
					trace[j].Center 		= 0;
					trace[j].Clutter 		= 0;
					trace[j].Num_in_Res 	= 0;
					trace[j].Res_i 			= 0;
				}

				// first mark:
				trace[0].Res_i 		= i;
				trace[0].Num_in_Res = k;
				trace[0].Center 	= Results[i-1].Centers[k-1];
				trace_length 		= 1;
				isClutter 			= 0;

				i2 = i;
				k2 = k;

				// try to continue the temp trace:
				while (1) {

					//if ((trace_length > 1) && (isClutter == 1)) 	isClutter_prev = isClutter;
					//else 											isClutter_prev = 0;

					i_prev = i2;
					k_prev = k2;


					tmp = GetNextClosest (&Results[i_prev-1], k_prev, &k2, &isClutter, deltas1, deltas2, clutterThr);
					i2 += tmp;

					if (tmp > 0) {

						if (Results[i2-1].isInTrace[k2-1] > 0) { 								// traces collision

							SetWholeTraceTo (Results, trace, trace_length, -1);
							break;

						} else if (trace_length >= Ncptrd) { 								// temp trace have more than Ncptrd points - append this to the array of traces

							if (*Ntraces >= MaxTracks) {
								fprintf(stderr, "Too many traces have been detected!\n");
								flag_too_many_traces = 1;
								break;
							} else {

								// append to the array:
								*Ntraces = *Ntraces + 1;
								for (j = 0; j < Nstep; j++) {
									memcpy (&traces[j][*Ntraces-1], &trace[j], sizeof(st_trace));
								}
								SetWholeTraceTo (Results, trace, trace_length, *Ntraces);

								// a search of new points for this array:
								while (1) {

									//isClutter_prev = isClutter;
									i_prev = i2;
									k_prev = k2;

									tmp = GetNextClosest (&Results[i_prev-1], k_prev, &k2, &isClutter, deltas1, deltas2, clutterThr);
									i2 += tmp;

									if (tmp > 0) {

										// find:
										for (j = 0; j < Nstep; j++) {
											if (traces[j][*Ntraces-1].Res_i == 0) break;
										}
										tmp = j;

										// append:
										traces[tmp][*Ntraces-1].Res_i 		= i2;
										traces[tmp][*Ntraces-1].Num_in_Res 	= k2;
										traces[tmp][*Ntraces-1].Center 		= Results[i2-1].Centers[k2-1];
										traces[tmp][*Ntraces-1].Clutter 	= isClutter;
										Results[i2-1].isInTrace[k2-1] 		= *Ntraces;

										// correction if clutter:
										if (isClutter) {
											traces[tmp][*Ntraces-1].Center = Soft_Part_II_interp1 (traces, tmp, *Ntraces-1, Nlin);
										}

										// also append point to temp trace:
										trace_length++;
										trace[trace_length].Res_i 		= i2;
										trace[trace_length].Num_in_Res 	= k2;
										trace[trace_length].Center 		= traces[tmp][*Ntraces-1].Center;
										trace[trace_length].Clutter 	= isClutter;

										// check end of data:
										if (i2 >= Nstep-3) {
											// ShowMeTheTrace(Results,traces,Ntraces) in MATLAB
											break;
										}

									} else {
										// ShowMeTheTrace(Results,traces,Ntraces) in MATLAB
										break;
									}

								}

								break;

							}

						} else { 															// append point to temp trace

							trace[trace_length].Res_i 		= i2;
							trace[trace_length].Num_in_Res 	= k2;
							trace[trace_length].Clutter 	= isClutter;
							trace[trace_length].Center 		= Results[i2-1].Centers[k2-1];
							trace_length++;

							if (i2 >= Nstep - 2) {
								SetWholeTraceTo (Results, trace, trace_length, -1);
								break;
							}

						}


					} else {
						SetWholeTraceTo (Results, trace, trace_length, -1);
						break;
					}

				}

				if (flag_too_many_traces) break;

			}

		}
		if (flag_too_many_traces) break;
	}


	// cleaning:
	free (trace);
	/*for (i = 0; i < Nstep; i++) free(traces[i]);
	free (traces);*/

	return traces;

}


float Soft_Part_II_interp1 (st_trace **traces, uint stepnum, uint tracenum, uint Nlin) {

	uint 	i;
	uint 	xi 		= traces[stepnum][tracenum].Res_i;
	uint 	*x 		= malloc(sizeof(uint)  * Nlin);
	float 	*y 		= malloc(sizeof(float) * Nlin);
	float 	*slope 	= malloc(sizeof(float) * (Nlin-1));
	float 	yi;

	for (i = 0; i < Nlin; i++) {
		x[i] = traces[stepnum - 3 + i][tracenum].Res_i;
		y[i] = traces[stepnum - 3 + i][tracenum].Center;
	}


	for (i = 0; i < Nlin-1; i++) {
		slope[i] = (y[i+1] - y[i]) / (x[i+1] - x[i]);
	}

	yi = (xi - x[Nlin-1]) * slope[Nlin-2] + y[Nlin-1];

	free (x);
	free (y);
	free (slope);
	return yi;

}



void Soft_Part_II_filtering (st_Results *Results, uint Nstep, float deltas1, float deltas2) {

	uint i, j, k;

	st_Results *R;
	st_Results *R1;
	st_Results *R2;

	for (i = 0; i < Nstep-2; i++) {

		// 3.1. neighbors frames:
		R 	= Results + i;
		R1 	= Results + i + 1;
		R2 	= Results + i + 2;

		for (j = 0; j < R->Num; j++) {

			// 3.2. check cross with any of R1:
			for (k = 0; k < R1->Num; k++) {
				if ((R->Begs[j] * (1 - deltas1) > R1->Fins[k]) || (R->Fins[j] * (1 + deltas1) < R1->Begs[k])) {
					continue;
				} else {
					R->isInTrace[j] = 0;
					break;
				}
			}
			if (R->isInTrace[j] == 0) continue;

			// 3.3. check cross with any of R1:
			for (k = 0; k < R2->Num; k++) {
				if ((R->Begs[j] * (1 - deltas2) > R2->Fins[k]) || (R->Fins[j] * (1 + deltas2) < R2->Begs[k])) {
					continue;
				} else {
					R->isInTrace[j] = 0;
					break;
				}
			}

		}

	}

}

uint GetNextClosest (st_Results *R, uint k, uint *n2, uchar *isClutter, float delta1, float delta2, float ClutterThr) {

	uint 	i2 			= 0;
	uint 	*found;
	uint 	found_i 	= 0;

	uint 	Begn 		= R[0].Begs[k-1];
	uint 	Finh 		= R[0].Fins[k-1];
	uint 	Wdth 		= R[0].Widths[k-1];
	uint 	k1;
	uint 	Begs1, Finh1;
	uint 	j;
	float 	ctmp, ctmp2;

	*isClutter 	= 0;
	*n2 		= 0;

	// For R+1:
	found = malloc (sizeof(uint) * R[1].Num);

	for (k1 = 0; k1 < R[1].Num; k1++) {
		if (R[1].isInTrace[k1] == 0) {
			Begs1 = R[1].Begs[k1];
			Finh1 = R[1].Fins[k1];
			if ( (((float) Begn * (1 - delta1)) > Finh1) || (((float) Finh * (1+delta1)) < Begs1) ) {

			} else {
				found[found_i++] = k1+1;
			}
		}
	}

	if (found_i > 0) {

		i2 = 1;

		if (found_i == 1) { 		// only one point - it's simple
			*n2 = found[0];
		} else {

			// min(abs(...)):
			ctmp 	= fabs(R[1].Centers[0] - R[0].Centers[k-1]);
			*n2 	= 1;
			for (j = 1; j < R[1].Num; j++) {
				ctmp2 = fabs(R[1].Centers[j] - R[0].Centers[k-1]);
				if (ctmp2 < ctmp) {
					*n2 = j+1;
				}
				ctmp = ctmp2;
			}

		}

		if (((float)R[1].Widths[*n2 - 1])/((float)Wdth) >= ClutterThr) {
			*isClutter = 1;
		}

		free (found);
		return i2;

	}
	free (found);


	// For R+2:
	found 	= malloc (sizeof(uint) * R[2].Num);
	found_i = 0;

	for (k1 = 0; k1 < R[2].Num; k1++) {
		if (R[2].isInTrace[k1] == 0) {
			Begs1 = R[2].Begs[k1];
			Finh1 = R[2].Fins[k1];
			if ( (((float) Begn * (1 - delta2)) > Finh1) || (((float) Finh * (1+delta2)) < Begs1) ) {

			} else {
				found[found_i++] = k1+1;
			}
		}
	}

	if (found_i > 0) {

		i2 = 2;

		if (found_i == 1) { 		// only one point - it's simple
			*n2 = found[0];
		} else {

			// min(abs(...)):
			ctmp 	= fabs(R[2].Centers[0] - R[0].Centers[k-1]);
			*n2 	= 1;
			for (j = 1; j < R[2].Num; j++) {
				ctmp2 = fabs(R[2].Centers[j] - R[0].Centers[k-1]);
				if (ctmp2 < ctmp) {
					*n2 = j+1;
				}
				ctmp = ctmp2;
			}

		}

		if (((float)R[2].Widths[*n2 - 1])/((float)Wdth) >= ClutterThr) {
			*isClutter = 1;
		}

		free (found);

		return i2;

	}



	free (found);

	return i2;

}

void SetWholeTraceTo (st_Results *R, st_trace *trace, const uint length, const int num) {

	uint i;

	for (i = 0; i < length; i++) {
		R[trace[i].Res_i - 1].isInTrace[trace[i].Num_in_Res - 1] = num;
	}

}

