#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "Soft_Part_I.h"
#include "Soft_Part_II.h"
#include "Soft_Part_III.h"
#include "polyfit.h"


void array_parameters (double *src, uint len, float *mean, float *min, float *max);
void polyval (double *Y, double *x, uint len, double *coef, uint degree);
double standard_deviation (double *src, uint len);


st_TracksInfo *Soft_Part_III (st_Results *Results, st_trace **traces, uint Nstep, uint Ntraces, uint *Lengths) {


	st_TracksInfo *TracksInfo = malloc (sizeof(st_TracksInfo) * Ntraces);
	uint n_trk;
	uint the_last;
	uint *x;
	double *y, *xf, *Y, *err;
	uint *tmp;
	uint *wtdhs;
	uint j;
	float ymean, ymin, ymax;
	double p[11];
	double errstd, errmean;
	double trh_brk;

	float sum_wtdhs = 0;

	for (n_trk = 0; n_trk < Ntraces; n_trk++) {

		// 1. calc length of track:
		for (the_last = 0; the_last < Nstep; the_last++) {
			if (traces[the_last][n_trk].Res_i == 0) break;
		}
		Lengths[n_trk] = the_last;

		// 2. get broadband of marks:
		x 		= malloc (sizeof(uint)   * the_last);
		xf 		= malloc (sizeof(double) * the_last);
		tmp 	= malloc (sizeof(uint)   * the_last);
		wtdhs 	= malloc (sizeof(uint)   * the_last);

		for (j = 0; j < the_last; j++) {
			x[j] 		= traces[j][n_trk].Res_i;
			xf[j] 		= (double) x[j];
			tmp[j] 		= traces[j][n_trk].Num_in_Res;
			wtdhs[j] 	= Results[x[j]].Types[tmp[j]];
		}

		// 3. checking for Stalker - once of the broadband radar:
		sum_wtdhs = 0;
		for (j = 0; j < the_last; j++) sum_wtdhs += wtdhs[j];
		if ( ((sum_wtdhs / the_last) / the_last) > 0.5 ) {
			fprintf(stdout, "===== Final Result ====================================\n");
			fprintf(stdout, "Most probable radar: %s\n", "Stalker 34G");
			fprintf(stdout, "With probability: %f%%\n", ((sum_wtdhs / the_last) / the_last));
			free (x);
			free (xf);
			free (tmp);
			free (wtdhs);
			free (TracksInfo);
			free (Lengths);
			return NULL;
		}

		// 4. get track values:
		y = malloc (sizeof(double) * the_last);
		Y = malloc (sizeof(double) * the_last);
		err = malloc (sizeof(double) * the_last);
		for (j = 0; j < the_last; j++) y[j] = traces[j][n_trk].Center;

		// 5. a measurement of a base parameters:
		array_parameters (y, the_last, &ymean, &ymin, &ymax);

		TracksInfo[n_trk].Fmean 	= ymean;
		TracksInfo[n_trk].dFmax 	= ymax - ymin;
		TracksInfo[n_trk].Length 	= x[the_last-1] - x[0];

		//polyfit (xf, y, the_last, 2, p);
		polyfit (the_last, 3, xf, y, p);
		TracksInfo[n_trk].Curvature = p[2] / p[1];

		//polyfit (xf, y, the_last, 10, p);
		polyfit (the_last, 11, xf, y, p);
		polyval (Y, xf, the_last, p, 10);

		for (j = 0; j < the_last; j++) err[j] = fabs(Y[j] - y[j]);
		errstd = standard_deviation (err, the_last);
		TracksInfo[n_trk].MSE = errstd;

		errmean = 0;
		for (j = 0; j < the_last; j++) errmean += err[j];
		errmean = errmean / ((double) the_last);

		trh_brk = errmean + 3.5 * errstd;

		TracksInfo[n_trk].isBreak = 0;
		for (j = 0; j < the_last; j++) {
			if (err[j] > trh_brk) {
				TracksInfo[n_trk].isBreak = 1;
				// [~,i_max] = max(err); in MATLAB, but i_max never used
				break;
			}
		}


		free (x);
		free (xf);
		free (tmp);
		free (wtdhs);
		free (y);
		free (Y);
		free (err);

	}

	return TracksInfo;

}

void array_parameters (double *src, uint len, float *mean, float *min, float *max) {

	uint i;

	*mean 	= 0;
	*min 	= src[0];
	*max 	= src[0];

	for (i = 0; i < len; i++) {
		*mean = *mean + src[i];
		if (src[i] < *min) *min = src[i];
		if (src[i] > *max) *max = src[i];
	}

	*mean = *mean / ((float) len);

}


void polyval (double *Y, double *x, uint len, double *coef, uint degree) {

	uint i, j;
	double tmp;

	for (i = 0; i < len; i++) {

		tmp = coef[0];
		for (j = 1; j <= degree; j++) {
			tmp += coef[j] * pow(x[i], j);
		}

		Y[i] = tmp;

	}


}


double standard_deviation (double *src, uint len) {

	uint i;
	double mean = 0;
	double diff = 0;
	double std;
	double N = (double) len;

	for (i = 0; i < len; i++) mean += src[i];
	mean = mean/N;

	for (i = 0; i < len; i++) diff += pow(fabs(src[i] - mean), 2);

	std = sqrt( (1/(N-1)) * diff);

	return std;

}


