#include <stdlib.h>
#include <stdio.h>
#include "stdafx.h"
#include "Soft_Part_I.h"
#include "Soft_Part_II.h"
#include "Soft_Part_III.h"
#include "recognition.h"

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>		/* open */
#include <unistd.h>		/* exit */
#include <time.h>
#include <sys/ioctl.h>	/* ioctl */
#include <sys/mman.h>
#include <string.h>



st_settings 	s;
st_fpgadata 	*fpgadata;
st_Results 		*Results;
st_trace 		**traces;
st_TracksInfo 	*TracksInfo;
uint 			*TracksLengths;

int get_Nstep 	(char *datafile);
int get_data 	(uint Nstep, char *datafile, st_fpgadata *fpgadata);

// int main (int argc, char *argv[]) {
int main()
{
	char *datafile;
	uint Nstep, Ntraces;
	uint i;
	char *Type;
	float Prob_of_recogn;

	// settings:
	s.Nwindow 		= FFT_WINDOW;
	s.Width 		= 40;
	s.delta1 		= 0.25;
	s.delta2 		= 0.5;
	s.Ncptrd 		= 30;
	s.MaxTracks 	= 30;
	s.MaxTraces 	= 5;
	s.clutterTrh 	= 3.0;
	s.Nlin 			= 3;


	// 1. get data from input file:
	Nstep = get_Nstep("fft_2_dma.bin");
	if (Nstep <= 0) {
		fprintf(stderr, "break input file\n");
		return -1;
	} else {
		fprintf(stdout, "number of Nstep: %d\n", Nstep);
	}


	fpgadata = malloc (sizeof(st_fpgadata) * Nstep);

	get_data(Nstep, datafile, fpgadata);
	printf("Got set of bins\n");

	// 2. Soft_part_I:
	Results = malloc(sizeof(st_Results) * Nstep);
	for (i = 0; i < Nstep; i++) {
		Soft_Part_I (&Results[i], fpgadata[i].Sf, fpgadata[i].Begs, fpgadata[i].Fins, s.Width);
	}

	// 3. Soft part II:
	traces = Soft_Part_II (Results, Nstep, &Ntraces, s.delta1, s.delta2, s.MaxTracks, s.Ncptrd, s.clutterTrh, s.Nlin);
	fprintf(stdout, "number of traces: %d\n", Ntraces);


	if (Ntraces > s.MaxTraces) {
		fprintf(stdout, "===== Final Result ====================================\n");
		fprintf(stdout, "Too many traces: unclassified.\n");
		free (fpgadata);
		free (Results);
		for (i = 0; i < Nstep; i++) free(traces[i]);
		free (traces);
		return -1;
	}

	// 6. Soft part III:
	 TracksLengths 	= malloc(sizeof(uint) * Ntraces);
	 TracksInfo 		= Soft_Part_III (Results, traces, Nstep, Ntraces, TracksLengths);

	// 8. Recognition:
	 Type = recognition (TracksInfo, TracksLengths, Ntraces, Nstep, traces, s.Nwindow, &Prob_of_recogn);

	// 9. Result:
	fprintf(stdout, "===== Final Result ====================================\n");
	fprintf (stdout, "Type: %s, Prob_of_recogn: %f\n", Type, Prob_of_recogn);

	// cleaning:
	free (fpgadata);
	free (Results);
	for (i = 0; i < Nstep; i++) free(traces[i]);
	free (traces);
	free (TracksInfo);
	free (TracksLengths);
	free (Type);


	return 0;

}
