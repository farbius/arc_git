#include <stdlib.h>
#include <string.h>
#include "Soft_Part_II.h"
#include "Soft_Part_III.h"
#include "recognition.h"

uchar IntersectEnough (st_trace **traces, uint *TracksLengths, uint Ntraces, uint Nstep);


char *recognition (st_TracksInfo *TracksInfo, uint *TracksLengths, uint Ntraces, uint Nstep, st_trace **traces, uint Nwindow, float *Prob_of_recogn) {

	char *Type = malloc(sizeof(char) * 32);
	float tlmax, tlmin;

	if (Ntraces == 1) {
		if (TracksInfo[0].dFmax < 10) {
			strcpy	(Type, "Unclassified");
			*Prob_of_recogn 	= 0.0;
		} else {
			if (TracksInfo[0].isBreak == 1) {
				strcpy	(Type, "Ramera_222");
				*Prob_of_recogn = TracksInfo[0].dFmax/(Nwindow/8);
			} else {
				if (TracksInfo[0].dFmax > 200) {
					strcpy	(Type, "Ramera_222");
					*Prob_of_recogn = 1-1/(TracksInfo[0].dFmax/25);
				} else {
					strcpy	(Type, "Unclassified");
					*Prob_of_recogn 	= 0.0;
				}
			}
		}
	} else {
		if (IntersectEnough(traces, TracksLengths, Ntraces, Nstep)) {
			strcpy	(Type, "Multa CD");
			if (TracksLengths[0] > TracksLengths[1]) {
				tlmax = TracksLengths[0];
				tlmin = TracksLengths[1];
			} else {
				tlmax = TracksLengths[1];
				tlmin = TracksLengths[0];
			}
			*Prob_of_recogn = tlmin / tlmax;
		} else {
			if ((Ntraces == 2) || (Ntraces == 3)) {
				strcpy	(Type, "Iskra");
				*Prob_of_recogn = ((float)(traces[0][1].Res_i) - (float)(traces[TracksLengths[0]-1][0].Res_i)) / ((float)Nstep);
			} else {
				strcpy	(Type, "Unclassified");
				*Prob_of_recogn 	= 0.0;
			}
		}

	}



	return Type;


}



uchar IntersectEnough (st_trace **traces, uint *TracksLengths, uint Ntraces, uint Nstep) {

	float LimLength 	= 0.9;
	//float LimIntersec 	= 0.5;

	uint i;
	//uint *Ilm = malloc(sizeof(uint) * Ntraces);
	uint longtracks = 0;

	for (i = 0; i < Ntraces; i++) {
		if (((float)TracksLengths[i])/((float)Nstep) >= LimLength) {
			//Ilm[longtracks] = i;
			longtracks++;
		}
	}

	if (longtracks < 2) {
		return 0;
	} else {
		// many calculation in MATLAB
		// ...

		return 1;
	}


}
